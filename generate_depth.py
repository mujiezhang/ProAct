#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import time
import shutil
from pathlib import Path


def run_command(cmd, shell=False, check=True):
    """Execute command; print original stderr/stdout to screen on failure"""
    # Do not capture output so child process writes directly to the terminal
    result = subprocess.run(cmd, shell=shell)
    if check and result.returncode != 0:
        # Exit with the same code without altering error messages
        sys.exit(result.returncode)
    return result



def check_tools():
    """Check if required tools are available"""
    required_tools = ['bbduk.sh', 'bwa', 'samtools']
    missing_tools = []
    
    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)
    
    if missing_tools:
        print("Error: The following tools are not found in PATH: {}".format(', '.join(missing_tools)))
        print("Ensure ProAct environment is activated.")
        sys.exit(1)



def quality_control(sample, read1, read2):
    """Quality control using bbduk.sh"""
    cmd = [
        'bbduk.sh',
        'in1={}'.format(read1),
        'in2={}'.format(read2),
        'out1={}_1_clean.fastq.gz'.format(sample),
        'out2={}_2_clean.fastq.gz'.format(sample),
        'ref=adapters,phix',
        'trimq=14',
        'maq=20',
        'maxns=1',
        'minlen=45'
    ]
    
    run_command(cmd)


def alignment(sample, ref_genome, threads):
    """Mapping using bwa"""
    # Check if index exists, create if not
    if not os.path.exists("{}.bwt".format(ref_genome)):
        index_cmd = ['bwa', 'index', ref_genome]
        run_command(index_cmd)
    
    # Mapping
    read1_clean = "{}_1_clean.fastq.gz".format(sample)
    read2_clean = "{}_2_clean.fastq.gz".format(sample)
    temp_bam = "{}_temp.bam".format(sample)
    
    # Use shell=True for pipe operations
    align_cmd = "bwa mem -a -t {} {} {} {} | samtools view -b -F 4 -o {}".format(
        threads, ref_genome, read1_clean, read2_clean, temp_bam
    )
    run_command(align_cmd, shell=True)


def filter_alignments(sample):
    """Filter alignments using pysam"""
    temp_bam = "{}_temp.bam".format(sample)
    filtered_bam = "{}_filtered.bam".format(sample)
    
    # Import pysam here to avoid import error if not installed
    try:
        import pysam
    except ImportError:
        print("Error: pysam is required for alignment filtering. Please install it with: pip install pysam")
        sys.exit(1)
    
    min_al, min_id, min_cov = 45, 0.97, 0.80
    
    with pysam.AlignmentFile(temp_bam, "rb") as in_bam:
        with pysam.AlignmentFile(filtered_bam, "wb", template=in_bam) as out_bam:
            for read in in_bam:
                if read.is_unmapped:
                    continue
                
                stats = read.get_cigar_stats()[0]
                al_len = sum(stats[:3])
                q_len = read.infer_read_length() or 0
                
                if al_len < min_al or q_len <= 0:
                    continue
                
                cov = (stats[0] + stats[1]) / q_len
                nm = read.get_tag("NM") if read.has_tag("NM") else al_len
                
                if (al_len - nm) / al_len >= min_id and cov >= min_cov:
                    out_bam.write(read)


def sort_and_depth(sample, threads):
    """Sort BAM and generate depth file"""
    filtered_bam = "{}_filtered.bam".format(sample)
    sorted_bam = "{}_sorted.bam".format(sample)
    depth_file = "{}.depth".format(sample)
    
    # Sort BAM
    sort_cmd = ['samtools', 'sort', '-@', str(threads), filtered_bam, '-o', sorted_bam]
    run_command(sort_cmd)
    
    # Generate depth
    depth_cmd = ['samtools', 'depth', '-aa', sorted_bam]
    with open(depth_file, 'w') as f:
        result = subprocess.run(depth_cmd, stdout=f, check=True)


def cleanup(sample, ref_genome):
    """Clean up intermediate files"""
    files_to_remove = [
        "{}_temp.bam".format(sample),
        "{}_filtered.bam".format(sample), 
        "{}_sorted.bam".format(sample),
        "{}_1_clean.fastq.gz".format(sample),
        "{}_2_clean.fastq.gz".format(sample)
    ]
    
    # BWA index files
    index_extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
    for ext in index_extensions:
        files_to_remove.append("{}.{}".format(ref_genome, ext))
    
    for file_path in files_to_remove:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
        except OSError:
            pass  # Ignore errors when removing files


def create_input_list(output_prefix, host_genome, raw_r1, raw_r2):
    """Derive sample prefix from read1 filename; do not create any file"""
    list_txt = None
    # Extract sample prefix from read1 filename
    sample_prefix = Path(raw_r1).stem
    # Remove common suffixes
    for suffix in ['_R1', '_1', '.fastq', '.fq']:
        if sample_prefix.endswith(suffix):
            sample_prefix = sample_prefix[:-len(suffix)]
            break
    return list_txt, sample_prefix


def process_sample(sample, ref_genome, read1, read2, threads):
    """Process a single sample through the complete pipeline"""
    # Stage 1: Quality Control
    quality_control(sample, read1, read2)
    
    # Stage 2: Alignment
    alignment(sample, ref_genome, threads)
    
    # Stage 3: Filter alignments
    filter_alignments(sample)
    
    # Stage 4: Sort and generate depth
    sort_and_depth(sample, threads)
    
    # Stage 5: Cleanup
    cleanup(sample, ref_genome)


def main():
    parser = argparse.ArgumentParser(
        description='Genome QC, mapping, filtering, sorting and depth generation (Python)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_depth.py -g host_genome.fna -1 raw_R1.fastq -2 raw_R2.fastq
  python generate_depth.py -g genome.fna -1 reads_R1.fq -2 reads_R2.fq -t 8 -o sample_001
        """
    )
    
    parser.add_argument('-g', '--genome', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('-1', '--read1', required=True,
                       help='Raw Read1 FASTQ')
    parser.add_argument('-2', '--read2', required=True,
                       help='Raw Read2 FASTQ')
    parser.add_argument('-t', '--threads', type=int, default=40,
                       help='Number of threads (default: 40)')
    parser.add_argument('-o', '--output-prefix', default='pipeline',
                       help='Output file prefix (default: pipeline)')
    
    args = parser.parse_args()
    
    # Check required files exist
    for file_path in [args.genome, args.read1, args.read2]:
        if not os.path.exists(file_path):
            print("Error: File not found: {}".format(file_path))
            sys.exit(1)
    
    # Check tools availability
    check_tools()
    # Create log files
    log_out = "{}.out".format(args.output_prefix)
    log_err = "{}.err".format(args.output_prefix)
    
    # Create input list and get sample name
    list_txt, sample_name = create_input_list(
        args.output_prefix, args.genome, args.read1, args.read2
    )
    
    start_time = time.time()

    process_sample(sample_name, args.genome, args.read1, args.read2, args.threads)
        
    end_time = time.time()
    elapsed = end_time - start_time

if __name__ == '__main__':
    main()
