#!/usr/bin/env python3
"""
Extract GTDB marker genes from a genome.

Provides both a functional API (run_extract_gtdb_mg) and a CLI entry point.
Adapted from GTDB-Tk identify module.
"""
import sys
import os
import argparse
import subprocess
from typing import Optional

from proact.extract_gtdb_mg.lib.Scan.PfamScan import PfamScan

sys.setrecursionlimit(2000)

hmmsearch_cmd = 'hmmsearch'
prodigal_cmd = 'prodigal'

# Marker information
BAC120_MARKERS = {
    "PFAM": [
        "PF00380.20.hmm", "PF00410.20.hmm", "PF00466.21.hmm",
        "PF01025.20.hmm", "PF02576.18.hmm", "PF03726.15.hmm",
    ],
    "TIGRFAM": [
        "TIGR00006.HMM", "TIGR00019.HMM", "TIGR00020.HMM",
        "TIGR00029.HMM", "TIGR00043.HMM", "TIGR00054.HMM",
        "TIGR00059.HMM", "TIGR00061.HMM", "TIGR00064.HMM",
        "TIGR00065.HMM", "TIGR00082.HMM", "TIGR00083.HMM",
        "TIGR00084.HMM", "TIGR00086.HMM", "TIGR00088.HMM",
        "TIGR00090.HMM", "TIGR00092.HMM", "TIGR00095.HMM",
        "TIGR00115.HMM", "TIGR00116.HMM", "TIGR00138.HMM",
        "TIGR00158.HMM", "TIGR00166.HMM", "TIGR00168.HMM",
        "TIGR00186.HMM", "TIGR00194.HMM", "TIGR00250.HMM",
        "TIGR00337.HMM", "TIGR00344.HMM", "TIGR00362.HMM",
        "TIGR00382.HMM", "TIGR00392.HMM", "TIGR00396.HMM",
        "TIGR00398.HMM", "TIGR00414.HMM", "TIGR00416.HMM",
        "TIGR00420.HMM", "TIGR00431.HMM", "TIGR00435.HMM",
        "TIGR00436.HMM", "TIGR00442.HMM", "TIGR00445.HMM",
        "TIGR00456.HMM", "TIGR00459.HMM", "TIGR00460.HMM",
        "TIGR00468.HMM", "TIGR00472.HMM", "TIGR00487.HMM",
        "TIGR00496.HMM", "TIGR00539.HMM", "TIGR00580.HMM",
        "TIGR00593.HMM", "TIGR00615.HMM", "TIGR00631.HMM",
        "TIGR00634.HMM", "TIGR00635.HMM", "TIGR00643.HMM",
        "TIGR00663.HMM", "TIGR00717.HMM", "TIGR00755.HMM",
        "TIGR00810.HMM", "TIGR00922.HMM", "TIGR00928.HMM",
        "TIGR00959.HMM", "TIGR00963.HMM", "TIGR00964.HMM",
        "TIGR00967.HMM", "TIGR01009.HMM", "TIGR01011.HMM",
        "TIGR01017.HMM", "TIGR01021.HMM", "TIGR01029.HMM",
        "TIGR01032.HMM", "TIGR01039.HMM", "TIGR01044.HMM",
        "TIGR01059.HMM", "TIGR01063.HMM", "TIGR01066.HMM",
        "TIGR01071.HMM", "TIGR01079.HMM", "TIGR01082.HMM",
        "TIGR01087.HMM", "TIGR01128.HMM", "TIGR01146.HMM",
        "TIGR01164.HMM", "TIGR01169.HMM", "TIGR01171.HMM",
        "TIGR01302.HMM", "TIGR01391.HMM", "TIGR01393.HMM",
        "TIGR01394.HMM", "TIGR01510.HMM", "TIGR01632.HMM",
        "TIGR01951.HMM", "TIGR01953.HMM", "TIGR02012.HMM",
        "TIGR02013.HMM", "TIGR02027.HMM", "TIGR02075.HMM",
        "TIGR02191.HMM", "TIGR02273.HMM", "TIGR02350.HMM",
        "TIGR02386.HMM", "TIGR02397.HMM", "TIGR02432.HMM",
        "TIGR02729.HMM", "TIGR03263.HMM", "TIGR03594.HMM",
        "TIGR03625.HMM", "TIGR03632.HMM", "TIGR03654.HMM",
        "TIGR03723.HMM", "TIGR03725.HMM", "TIGR03953.HMM",
    ],
}

AR53_MARKERS = {
    "PFAM": [
        "PF04919.13.hmm", "PF07541.13.hmm", "PF01000.27.hmm",
        "PF00687.22.hmm", "PF00466.21.hmm", "PF00827.18.hmm",
        "PF01280.21.hmm", "PF01090.20.hmm",
        "PF01200.19.hmm", "PF01015.19.hmm", "PF00900.21.hmm", "PF00410.20.hmm",
    ],
    "TIGRFAM": [
        "TIGR00037.HMM", "TIGR00064.HMM", "TIGR00111.HMM",
        "TIGR00134.HMM", "TIGR00279.HMM", "TIGR00291.HMM", "TIGR00323.HMM",
        "TIGR00335.HMM", "TIGR00373.HMM", "TIGR00405.HMM", "TIGR00448.HMM",
        "TIGR00483.HMM", "TIGR00491.HMM", "TIGR00522.HMM", "TIGR00967.HMM",
        "TIGR00982.HMM", "TIGR01008.HMM", "TIGR01012.HMM", "TIGR01018.HMM",
        "TIGR01020.HMM", "TIGR01028.HMM", "TIGR01046.HMM", "TIGR01052.HMM",
        "TIGR01171.HMM", "TIGR01213.HMM", "TIGR01952.HMM", "TIGR02236.HMM",
        "TIGR02338.HMM", "TIGR02389.HMM", "TIGR02390.HMM", "TIGR03626.HMM",
        "TIGR03627.HMM", "TIGR03628.HMM", "TIGR03629.HMM", "TIGR03670.HMM",
        "TIGR03671.HMM", "TIGR03672.HMM", "TIGR03673.HMM", "TIGR03674.HMM",
        "TIGR03676.HMM", "TIGR03680.HMM",
    ],
}

marker_dict = dict()
for marker in BAC120_MARKERS['PFAM'] + BAC120_MARKERS['TIGRFAM']:
    marker_dict[marker[:-4]] = 'BAC'
for marker in AR53_MARKERS['PFAM'] + AR53_MARKERS['TIGRFAM']:
    if marker in marker_dict:
        marker_dict[marker[:-4]] = 'BOTH'
    else:
        marker_dict[marker[:-4]] = 'AR'


def has_tool(cmd):
    """Check if a command-line tool is available in PATH."""
    try:
        env = os.environ.copy()
        proc = subprocess.Popen(['which', cmd], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, env=env, encoding='utf-8')
        output, _ = proc.communicate()
        return len(output) > 0
    except Exception:
        return False


# --- Run commands ---

def run_prodigal(genome_file, tmp_folder, closed_ends=False, meta=False):
    closed_ends_flag = '-c' if closed_ends else ''
    genome_basename = os.path.basename(genome_file)
    genome_basename = genome_basename.rsplit('.', 1)[0]

    aa_out = os.path.join(tmp_folder, "{}.faa".format(genome_basename))
    nt_out = os.path.join(tmp_folder, "{}.fna".format(genome_basename))

    procedure = 'meta' if meta else 'single'
    args = ['prodigal', '-m', '-p', procedure, '-q', '-g', '11',
            '-a', aa_out, '-d', nt_out, '-i', genome_file]

    print(' '.join(args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        print('Non-zero exit code returned when running prodigal: {stdout}')
        return False, None, None

    return True, aa_out, nt_out


def run_pfam_search(genes_aa, hmm_dir, tmp_folder, cpu=1):
    genes_basename = os.path.basename(genes_aa)
    genes_basename = genes_basename.rsplit('.', 1)[0]

    pfam_output = os.path.join(tmp_folder, "{}_pfam.hits".format(genes_basename))

    pfam_scan = PfamScan(cpu=cpu, fasta=genes_aa, dir=hmm_dir)
    pfam_scan.search()
    pfam_scan.write_results(pfam_output, None, None, None, None)

    return True, pfam_output


def run_tigrfam_search(genes_aa, hmm_file, tmp_folder, cpu=1):
    genes_basename = os.path.basename(genes_aa)
    genes_basename = genes_basename.rsplit('.', 1)[0]

    tigr_output = os.path.join(tmp_folder, "{}_tigr".format(genes_basename))
    tigr_output_hits = os.path.join(tmp_folder, "{}_tigr.hits".format(genes_basename))

    args = [hmmsearch_cmd, '-o', tigr_output, '--tblout', tigr_output_hits,
            '--noali', '--notextw', '--cut_nc', '--cpu', str(cpu), hmm_file, genes_aa]

    p = subprocess.Popen(args, stdout=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        print('Non-zero exit code returned when running prodigal: {stdout}')
        return False, None, None

    return True, tigr_output_hits


# --- Extract commands ---

class Hit(object):
    """HMM hit with gene_id, hmm_id, e_val, bit_score."""

    def __init__(self, gene_id, hmm_id, e_val, bit_score):
        self.gene_id = gene_id
        self.hmm_id = hmm_id
        self.e_val = e_val
        self.bit_score = bit_score

    def __repr__(self):
        return f'{self.gene_id} {self.hmm_id} ({self.e_val}/{self.bit_score})'

    def __eq__(self, other):
        return (isinstance(other, Hit) and other.gene_id == self.gene_id
                and other.hmm_id == self.hmm_id and other.e_val == self.e_val
                and other.bit_score == self.bit_score)

    def __lt__(self, other):
        if self.bit_score < other.bit_score:
            return True
        elif self.bit_score == other.bit_score:
            if self.e_val > other.e_val:
                return True
            elif self.e_val == other.e_val:
                if self.hmm_id > other.hmm_id:
                    return True
                elif self.hmm_id == other.hmm_id:
                    return self.gene_id > other.gene_id
        return False

    def __gt__(self, other):
        raise NotImplemented

    def __hash__(self):
        return hash(f'{self.gene_id}_{self.hmm_id}_{self.e_val}_{self.bit_score}')

    def hmm_str(self):
        return f'{self.hmm_id},{self.e_val},{self.bit_score}'


def add_hit(hit_dict, gene_id, hmm_id, e_val, bit_score):
    if gene_id not in hit_dict:
        hit_dict[gene_id] = dict()
    new_hit = Hit(gene_id, hmm_id, e_val, bit_score)
    if hmm_id in hit_dict[gene_id]:
        if hit_dict[gene_id][hmm_id] < new_hit:
            hit_dict[gene_id][hmm_id] = new_hit
    else:
        hit_dict[gene_id][hmm_id] = new_hit


def get_top_hit(hit_dict, gene_id):
    if gene_id not in hit_dict or len(hit_dict[gene_id]) == 0:
        return None
    return sorted(hit_dict[gene_id].values(), reverse=True)[0]


def contains_gene_id(hit_dict, gene_id):
    return gene_id in hit_dict


def add_hit_tigrfam(hit_dict, gene_id, hmm_id, e_val, bit_score):
    if contains_gene_id(hit_dict, gene_id):
        if get_top_hit(hit_dict, gene_id) < Hit(gene_id, hmm_id, e_val, bit_score):
            hit_dict[gene_id] = dict()
            add_hit(hit_dict, gene_id, hmm_id, e_val, bit_score)
    else:
        add_hit(hit_dict, gene_id, hmm_id, e_val, bit_score)


def write_tophits(hit_dict, output):
    header = ['Gene Id', 'Top hits (Family id,e-value,bitscore)']
    with open(output, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')
        for gene_id, hits in sorted(hit_dict.items()):
            out_hits = []
            for cur_hit in sorted(hits.values(), reverse=True):
                out_hits.append(cur_hit.hmm_str())
            concat_hits = ';'.join(out_hits)
            outfile.write(f'{gene_id}\t{concat_hits}\n')


def tophit_pfam(pfam_file):
    hit_dict = dict()
    with open(pfam_file, 'r') as fh_pfam:
        for line in fh_pfam:
            if line[0] == '#' or not line.strip():
                continue
            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            add_hit(hit_dict, gene_id, hmm_id, evalue, bitscore)
    return hit_dict


def tophit_tigr(tigrfam_file):
    hit_dict = dict()
    with open(tigrfam_file, 'r') as fh_tigrfam:
        for line in fh_tigrfam:
            if line[0] == '#':
                continue
            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[3]
            evalue = float(line_split[4])
            bitscore = float(line_split[5])
            add_hit_tigrfam(hit_dict, gene_id, hmm_id, evalue, bitscore)
    return hit_dict


class Record:
    def __init__(self, header):
        self.header = header
        self.sequence = ""

    def append(self, sequence):
        self.sequence += sequence


def ReadFasta(file_path):
    header = ""
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):
                header = line
                if sequence and header:
                    record = Record(header)
                    record.sequence = sequence
                    yield record
                sequence = ""
            else:
                sequence += line
    if sequence and header:
        record = Record(header)
        record.sequence = sequence
        yield record
    return None


def get_header_to_mg(top_hits_pfam_filepath, top_hits_tigr_filepath):
    hit_dictionary = dict()
    for file_path in (top_hits_pfam_filepath, top_hits_tigr_filepath):
        with open(file_path, 'r') as file:
            next(file, None)
            for line in file:
                line = line.rstrip()
                if not line:
                    continue
                tokens = line.split('\t')
                gene_id = tokens[0]
                try:
                    hmm_id, _, _ = tokens[1].split(',')
                except (IndexError, ValueError):
                    continue
                if hmm_id not in hit_dictionary:
                    hit_dictionary[hmm_id] = list()
                hit_dictionary[hmm_id].append(gene_id)

    extract_list = dict()
    for marker_name, gene_id_list in hit_dictionary.items():
        for gene_id in gene_id_list:
            extract_list[gene_id] = marker_name
    return extract_list


def write_markers(output_path, header2mg, sequence_file):
    with open(output_path, 'w') as output:
        for record in ReadFasta(sequence_file):
            gene_name = record.header[1:].split(' # ')[0]
            if gene_name in header2mg:
                record.header = ">" + header2mg[gene_name]
                output.write(record.header + '\n')
                output.write(record.sequence + '\n')


def make_dir(path, silent=True):
    if os.path.isdir(path):
        return
    try:
        os.makedirs(path)
    except OSError:
        if not silent:
            print("Creation of the directory %s failed" % path)


def rm_dir(path, silent=True):
    try:
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
        os.rmdir(path)
    except OSError:
        if not silent:
            print("Deletion of the directory %s failed" % path)


def _get_data_dir():
    """Return the path to the HMM data directory within the installed package."""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'hmm')


def run_extract_gtdb_mg(genome, output_dir, tmp_dir, threads=1,
                         overwrite_tmp=True, list_only=True, silent=True,
                         meta=False, keep_tmp=False):
    """
    Functional API to extract GTDB marker genes.

    Parameters
    ----------
    genome : str
        Genome file path.
    output_dir : str
        Marker gene output folder.
    tmp_dir : str
        Temporary folder.
    threads : int
        Number of threads.
    overwrite_tmp : bool
        Overwrite existing tmp directory.
    list_only : bool
        Output a tab-delimited file instead of sequences.
    silent : bool
        Suppress output.
    meta : bool
        Use prodigal meta mode.
    keep_tmp : bool
        Keep temporary files.
    """
    # Resolve HMM data paths from the installed package
    PFAM_FOLDER = os.path.join(_get_data_dir(), 'pfam')
    TIGRFAM_FILE = os.path.join(_get_data_dir(), 'tigrfam', 'tigrfam.hmm')

    if not has_tool(hmmsearch_cmd):
        print("Need tool {}".format(hmmsearch_cmd))
        sys.exit(1)
    if not has_tool(prodigal_cmd):
        print("Need tool {}".format(prodigal_cmd))
        sys.exit(1)

    if os.path.isdir(tmp_dir) and not overwrite_tmp:
        print("Tmp dir {} already exists. Please delete dir first".format(tmp_dir))
        sys.exit(1)

    make_dir(tmp_dir)

    marker_genes_file = "marker_genes.tsv"
    list_output = None

    if list_only:
        list_output = open(os.path.join(output_dir, marker_genes_file), 'w')
        list_output.write("gene_id\tmarker\tstart\tend\n")

    input_genomes = genome.split(',')
    input_len = len(input_genomes)

    for index in range(input_len):
        input_genome = input_genomes[index]
        input_basename = os.path.basename(input_genome[index])
        input_basename = input_basename.rsplit('.', 1)[0]

        success, aa_file, nt_file = run_prodigal(input_genome, tmp_dir, meta=meta)

        pfam_success, pfam_hits = run_pfam_search(aa_file, PFAM_FOLDER, tmp_dir, cpu=threads)
        tigr_success, tigr_hits = run_tigrfam_search(aa_file, TIGRFAM_FILE, tmp_dir, cpu=threads)

        dict_tigr = tophit_tigr(tigr_hits)
        dict_pfam = tophit_pfam(pfam_hits)

        tigr_tophit_file = os.path.join(tmp_dir, '{}_tigr_tophit.tsv'.format(input_basename))
        pfam_tophit_file = os.path.join(tmp_dir, '{}_pfam_tophit.tsv'.format(input_basename))

        write_tophits(dict_tigr, tigr_tophit_file)
        write_tophits(dict_pfam, pfam_tophit_file)

        header2mg = get_header_to_mg(pfam_tophit_file, tigr_tophit_file)

        mg_basename = os.path.join(output_dir, input_basename)

        if list_only:
            gene_id_to_header = {}
            with open(aa_file, 'r') as f_aa:
                for line in f_aa:
                    if line.startswith('>'):
                        full_header = line.strip()[1:]
                        gene_id = full_header.split(' # ')[0]
                        gene_id_to_header[gene_id] = full_header
            for header, mg in header2mg.items():
                full_header = gene_id_to_header.get(header)
                if full_header:
                    parts = full_header.split(' # ')
                    start_pos = parts[1]
                    end_pos = parts[2]
                    list_output.write("{}\t{}\t{}\t{}\n".format(header, mg, start_pos, end_pos))
                else:
                    list_output.write("{}\t{}\tNA\tNA\n".format(header, mg))
        else:
            write_markers("{}.faa".format(mg_basename), header2mg, aa_file)
            if nt_file:
                write_markers("{}.fna".format(mg_basename), header2mg, nt_file)

    if list_output:
        list_output.close()

    if not keep_tmp:
        rm_dir(tmp_dir)


def main():
    parser = argparse.ArgumentParser(
        prog='proact-extract-gtdb-mg',
        description='Extract GTDB marker genes from a genome',
    )

    parser.add_argument("-a", "--aa", type=str,
                        help="One or more files with genes (Protein)")
    parser.add_argument("-n", "--nt", type=str,
                        help="One or more files with genes (DNA)")
    parser.add_argument("-g", "--genome", type=str,
                        help="Genome file(s)")

    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="Marker_gene output folder")
    parser.add_argument("-p", "--tmp", type=str, required=True,
                        help="Temporary folder")
    parser.add_argument("-k", "--keep_tmp", action="store_true",
                        help="Keep temporary files")
    parser.add_argument("-x", "--overwrite_tmp", action="store_true",
                        help="Ignore if tmp dir already exists")
    parser.add_argument("-t", "--threads", type=int,
                        help="Number of threads to use.")
    parser.add_argument("-l", "--list_only", action="store_true",
                        help="Do not output marker genes as sequence. "
                             "Output a single tab-delimited file.")
    parser.add_argument("-s", "--silent", action="store_true",
                        help="No output")
    parser.add_argument("-m", "--meta", action="store_true",
                        help="Input is metagenome (proteins from multiple potential genomes)")

    args = parser.parse_args()

    from_genomes = True

    if (args.aa or args.nt) and args.genome:
        print("options (--aa|--nt) and --genomes are mutually exclusive.")
        sys.exit(1)

    if not args.genome and not args.aa:
        print("Either --genome or --aa must be specified")
        sys.exit(1)

    input_len = 0
    if args.aa:
        from_genomes = False
        input_aa = args.aa.split(',')
        input_nt = None
        if args.nt:
            input_nt = args.nt.split(',')
            if len(input_aa) != len(input_nt):
                print("If --aa and --nt are specified, both must be of same length")
                sys.exit(1)
        input_len = len(input_aa)

    if args.genome:
        from_genomes = True
        input_genomes = args.genome.split(',')
        input_len = len(input_genomes)

    threads = args.threads if args.threads else 1

    # Resolve HMM data paths from the installed package
    script_dir = os.path.dirname(os.path.realpath(__file__))
    PFAM_FOLDER = os.path.join(script_dir, 'hmm', 'pfam')
    TIGRFAM_FILE = os.path.join(script_dir, 'hmm', 'tigrfam', 'tigrfam.hmm')

    if not has_tool(hmmsearch_cmd):
        print("Need tool {}".format(hmmsearch_cmd))
        sys.exit(1)
    if not has_tool(prodigal_cmd):
        print("Need tool {}".format(prodigal_cmd))
        sys.exit(1)

    output_dir = args.output_dir
    tmp_dir = args.tmp
    list_only = args.list_only
    silent = args.silent
    meta = args.meta

    marker_genes_file = "marker_genes.tsv"
    list_output = None

    if os.path.isdir(tmp_dir) and not args.overwrite_tmp:
        print("Tmp dir {} already exists. Please delete dir first".format(tmp_dir))
        sys.exit(1)

    make_dir(tmp_dir)

    if list_only:
        list_output = open(os.path.join(output_dir, marker_genes_file), 'w')
        list_output.write("gene_id\tmarker\tstart\tend\n")

    for index in range(input_len):
        aa_file = None
        nt_file = None

        if from_genomes:
            input_genome = input_genomes[index]
            input_basename = os.path.basename(input_genome)
            input_basename = input_basename.rsplit('.', 1)[0]
            success, aa_file, nt_file = run_prodigal(input_genome, tmp_dir, meta=meta)
        else:
            aa_file = input_aa[index]
            nt_file = input_nt[index] if input_nt else None
            input_basename = os.path.basename(aa_file)
            input_basename = input_basename.rsplit('.', 1)[0]

        pfam_success, pfam_hits = run_pfam_search(aa_file, PFAM_FOLDER, tmp_dir, cpu=threads)
        tigr_success, tigr_hits = run_tigrfam_search(aa_file, TIGRFAM_FILE, tmp_dir, cpu=threads)

        dict_tigr = tophit_tigr(tigr_hits)
        dict_pfam = tophit_pfam(pfam_hits)

        tigr_tophit_file = os.path.join(tmp_dir, '{}_tigr_tophit.tsv'.format(input_basename))
        pfam_tophit_file = os.path.join(tmp_dir, '{}_pfam_tophit.tsv'.format(input_basename))

        write_tophits(dict_tigr, tigr_tophit_file)
        write_tophits(dict_pfam, pfam_tophit_file)

        header2mg = get_header_to_mg(pfam_tophit_file, tigr_tophit_file)

        mg_basename = os.path.join(output_dir, input_basename)

        if list_only:
            gene_id_to_header = {}
            with open(aa_file, 'r') as f_aa:
                for line in f_aa:
                    if line.startswith('>'):
                        full_header = line.strip()[1:]
                        gene_id = full_header.split(' # ')[0]
                        gene_id_to_header[gene_id] = full_header
            for header, mg in header2mg.items():
                full_header = gene_id_to_header.get(header)
                if full_header:
                    parts = full_header.split(' # ')
                    start_pos = parts[1]
                    end_pos = parts[2]
                    list_output.write("{}\t{}\t{}\t{}\n".format(header, mg, start_pos, end_pos))
                else:
                    list_output.write("{}\t{}\tNA\tNA\n".format(header, mg))
        else:
            write_markers("{}.faa".format(mg_basename), header2mg, aa_file)
            if nt_file:
                write_markers("{}.fna".format(mg_basename), header2mg, nt_file)

    if list_output:
        list_output.close()

    if not args.keep_tmp:
        rm_dir(tmp_dir)


if __name__ == '__main__':
    main()
