#!/usr/bin/env python3
"""
ProAct Pipeline - Phage activity analysis

Description:
    End-to-end pipeline with four steps:
    1. QC, mapping, filtering, sorting, and depth generation
    2. Extract GTDB marker genes
    3. Compute depth statistics for marker genes and phage regions
    4. Compute PtoH

Requirements:
    - Python >=3.6
    - bwa, samtools, bbduk.sh
    - pysam (for alignment filtering)
Usage:
    python proact_pipeline.py \\
        -g <host_genome.fna> \\
        -1 <read1.fastq> \\
        -2 <read2.fastq> \\
        -p <phage_info.tsv> \\
        [-o <output_dir>] \\
        [-t <threads>] \\
        [--keep-tmp]

Version: 1.0
"""

import argparse
import os
import sys
import subprocess
import time
import shutil
from pathlib import Path


class ProActPipeline:
    def __init__(self, args):
        self.args = args
        self.start_time = time.time()
        
        # Set up directories
        self.output_dir = Path(args.output_dir)
        self.tmp_dir = self.output_dir / "tmp"
        
        # Create directories
        self.output_dir.mkdir(exist_ok=True)
        
        # Extract sample name from read1 filename
        self.sample_name = Path(args.read1).stem
        for suffix in ['_R1', '_1', '.fastq', '.fq']:
            if self.sample_name.endswith(suffix):
                self.sample_name = self.sample_name[:-len(suffix)]
                break
        
        print("="*60)
        print("ProAct Pipeline - Phage activity analysis")
        print("="*60)
        print("Sample: {}".format(self.sample_name))
        print("Output directory: {}".format(self.output_dir))
        print("Threads: {}".format(args.threads))
        print("Prodigal mode: {}".format("meta" if args.meta else "single"))
        print("="*60)


    def run_command(self, cmd, description="", check=True):
        """Execute command, suppress normal output; on failure print raw stdout/stderr"""
        if isinstance(cmd, str):
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True)
        if check and result.returncode != 0:
            # Print original outputs as-is, then exit with the same return code
            if result.stdout:
                print(result.stdout, end="")
            if result.stderr:
                print(result.stderr, end="")
            sys.exit(result.returncode)
        return result


    def check_dependencies(self):
        """Check if required tools are available"""
        required_tools = ['bwa', 'samtools', 'bbduk.sh', 'python']
        missing_tools = []
        
        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            print("ERROR: The following required tools are not found in PATH:")
            for tool in missing_tools:
                print("  - {}".format(tool))
            print("Please ensure ProAct environment is activated.")
            sys.exit(1)
        
        # Check Python modules
        try:
            import pysam
        except ImportError:
            print("ERROR: pysam module is required. Install with: pip install pysam")
            sys.exit(1)


    def step1_process_reads(self):
        """Step 1: Quality control, mapping, filtering, and depth calculation"""
        step_start = time.time()
        print("\n[STEP 1] QC, mapping, filtering, depth calculation...")
        
        depth_script = Path(__file__).parent / "generate_depth.py"
        if not depth_script.exists():
            print("ERROR: generate_depth.py not found")
            sys.exit(1)
        
        # Use absolute paths to avoid issues after chdir
        genome_abs = os.path.abspath(self.args.genome)
        read1_abs = os.path.abspath(self.args.read1)
        read2_abs = os.path.abspath(self.args.read2)
        cmd = [
            'python', str(depth_script),
            '-g', genome_abs,
            '-1', read1_abs,
            '-2', read2_abs,
            '-t', str(self.args.threads),
            '-o', self.sample_name
        ]
        
        # Change to output directory for this step
        original_cwd = os.getcwd()
        os.chdir(self.output_dir)
        
        try:
            self.run_command(cmd, "Quality control and mapping")
            
            # Check if depth file was created
            depth_file = Path("{}.depth".format(self.sample_name))
            if not depth_file.exists():
                print("ERROR: Depth file not created: {}".format(depth_file))
                sys.exit(1)
            
            step_time = time.time() - step_start
            print("[STEP 1] Done ({:.1f}s)".format(step_time))
            
        finally:
            os.chdir(original_cwd)


    def step2_extract_markers(self):
        """Step 2: Extract GTDB marker genes"""
        step_start = time.time()
        print("\n[STEP 2] Extract GTDB marker genes...")
        
        extract_script = Path(__file__).parent / "extract_gtdb_mg" / "extract_gtdb_mg.py"
        if not extract_script.exists():
            print("ERROR: extract_gtdb_mg.py not found")
            sys.exit(1)
        
        # Create tmp directory if it doesn't exist
        self.tmp_dir.mkdir(exist_ok=True)
        
        cmd = [
            'python', str(extract_script),
            '-g', self.args.genome,
            '-o', str(self.output_dir),
            '-p', str(self.tmp_dir),
            '-t', str(self.args.threads),
            '-x',
            '-l',  # list only mode
            '-s'   # silent mode
        ]
        
        # Add meta flag if specified
        if self.args.meta:
            cmd.append('-m')
        
        self.run_command(cmd, "GTDB marker gene extraction")
        
        # Check if marker genes file was created
        marker_file = self.output_dir / "marker_genes.tsv"
        if not marker_file.exists():
            print("ERROR: Marker genes file not created: {}".format(marker_file))
            sys.exit(1)
        
        step_time = time.time() - step_start
        print("[STEP 2] Done ({:.1f}s)".format(step_time))


    def step3_calculate_counts(self):
        """Step 3: Calculate marker gene and phage depth statistics"""
        step_start = time.time()
        print("\n[STEP 3] Calculate marker gene and phage depth statistics...")
        
        counts_script = Path(__file__).parent / "get_MG_and_phage_depth.py"
        if not counts_script.exists():
            print("ERROR: get_MG_and_phage_depth.py not found")
            sys.exit(1)
        
        depth_file = self.output_dir / "{}.depth".format(self.sample_name)
        marker_file = self.output_dir / "marker_genes.tsv"
        
        cmd = [
            'python', str(counts_script),
            self.args.phage_info,
            self.args.genome,
            str(depth_file),
            str(marker_file),
            str(self.output_dir)
        ]
        
        self.run_command(cmd, "Marker gene and phage depth calculation")
        
        # Check if output files were created
        expected_files = [
            self.output_dir / "marker_gene_counts.tsv",
            self.output_dir / "phage_counts.tsv",
            self.output_dir / "host_counts.tsv"
        ]
        
        for file_path in expected_files:
            if not file_path.exists():
                print("ERROR: Expected output file not created: {}".format(file_path))
                sys.exit(1)
        
        step_time = time.time() - step_start
        print("[STEP 3] Done ({:.1f}s)".format(step_time))


    def step4_calculate_ptoh(self):
        """Step 4: Calculate PtoH"""
        step_start = time.time()
        print("\n[STEP 4] Compute PtoH...")
        
        ptoh_script = Path(__file__).parent / "calculate_PtoH.py"
        if not ptoh_script.exists():
            print("ERROR: calculate_PtoH.py not found")
            sys.exit(1)
        
        phage_counts = self.output_dir / "phage_counts.tsv"
        host_counts = self.output_dir / "host_counts.tsv"
        ptoh_output = self.output_dir / "PtoH_results.tsv"
        
        cmd = [
            'python', str(ptoh_script),
            str(phage_counts),
            str(host_counts),
            str(ptoh_output)
        ]
        
        self.run_command(cmd, "PtoH calculation")
        
        if not ptoh_output.exists():
            print("ERROR: PtoH results file not created: {}".format(ptoh_output))
            sys.exit(1)
        
        step_time = time.time() - step_start
        print("[STEP 4] Done ({:.1f}s)".format(step_time))


    def cleanup(self):
        """Clean up temporary files if requested"""
        depth_file = self.output_dir / "{}.depth".format(self.sample_name)
        if depth_file.exists():
            try:
                depth_file.unlink()
            except OSError:
                pass
        
        if not self.args.keep_tmp:
            if self.tmp_dir.exists():
                shutil.rmtree(self.tmp_dir)
            
            # Remove intermediate files from output directory
            cleanup_patterns = ['*_temp.bam', '*_filtered.bam', '*_sorted.bam', 
                              '*_clean.fastq.gz', '*.bwt', '*.amb', '*.ann', '*.pac', '*.sa']
            
            for pattern in cleanup_patterns:
                for file_path in self.output_dir.glob(pattern):
                    try:
                        file_path.unlink()
                    except OSError:
                        pass


    def generate_summary(self):
        """Silent summary (no-op)"""
        return


    def run(self):
        """Run the complete pipeline"""
        try:
            self.check_dependencies()
            self.step1_process_reads()
            self.step2_extract_markers()
            self.step3_calculate_counts()
            self.step4_calculate_ptoh()
            self.cleanup()
            # Silent mode: no summary
            
        except KeyboardInterrupt:
            print("\n[WARNING] Pipeline interrupted by user")
            sys.exit(1)
        except Exception as e:
            print("\n[ERROR] Pipeline failed: {}".format(e))
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='ProAct Pipeline - Phage activity analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python proact_pipeline.py -g host_genome.fna -1 sample_R1.fastq -2 sample_R2.fastq -p phage_info.tsv

  # Specify output dir, threads and meta mode
  python proact_pipeline.py  -g genome.fna  -1 reads_R1.fq  -2 reads_R2.fq -p phage_regions.tsv -o results_sample1 -t 40 -m --keep-tmp

Outputs:
  marker_genes.tsv              - GTDB marker genes list
  marker_gene_counts.tsv        - depth stats for marker genes
  phage_counts.tsv              - depth stats for phage regions
  host_counts.tsv               - median depth of marker genes
  PtoH_results.tsv              - final prediction results
        """
    )
    
    # Required arguments
    parser.add_argument('-g', '--genome', required=True,
                       help='Host genome FASTA file')
    parser.add_argument('-1', '--read1', required=True,
                       help='Read1 FASTQ file')
    parser.add_argument('-2', '--read2', required=True,
                       help='Read2 FASTQ file')
    parser.add_argument('-p', '--phage-info', required=True,
                       help='A tab-delimited file with columns: virla_name, host_contig, start, end')
    
    # Optional arguments
    parser.add_argument('-o', '--output-dir', default='proact_results',
                       help='Output directory (default: proact_results)')
    parser.add_argument('-t', '--threads', type=int, default=40,
                       help='Number of threads (default: 40)')
    parser.add_argument('-m', '--meta', action='store_true',
                       help='Run prodigal in meta mode (for metagenomes)')
    parser.add_argument('--keep-tmp', action='store_true',
                       help='Keep temporary files')
    
    args = parser.parse_args()
    
    # Validate input files
    required_files = [args.genome, args.read1, args.read2, args.phage_info]
    for file_path in required_files:
        if not os.path.exists(file_path):
            print("ERROR: Input file not found: {}".format(file_path))
            sys.exit(1)
    
    # Run pipeline
    pipeline = ProActPipeline(args)
    pipeline.run()


if __name__ == '__main__':
    main()
