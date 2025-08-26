#!/usr/bin/env python3
"""
Script: get_MG_and_phage_depth.py
Description:
  Compute depth statistics for marker genes and phage regions of a single sample.

Usage:
  python get_MG_and_phage_depth.py \
    phage_info.tsv \
    host_genome.fna \
    sample.depth \
    marker_genes.tsv \
    ./result_dir

phage_info.tsv format:
  <viral_name>\t<host_contig>\t<start>\t<end>
"""
import os
import sys
import statistics

if len(sys.argv) != 6:
    print(__doc__)
    sys.exit(1)

phage_info_path, host_genome_path, depth_path, gene_path, output_dir = sys.argv[1:]

# Load phage_info (four columns: name, contig, start, end)
def load_phage_info(path):
    phage_list = []
    with open(path, 'r') as f:
        for line in f:
            raw = line.strip()
            if not raw:
                continue
            values = raw.split('\t')
            if len(values) < 4:
                continue
            name = values[0]
            contig = values[1]
            # Skip potential header
            try:
                start = int(values[2])
                end = int(values[3])
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            phage_list.append({
                'phage_name': name,
                'contig': contig,
                'start': start,
                'end': end
            })
    return phage_list

phage_list = load_phage_info(phage_info_path)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def load_depth(path):
    depth_list = []
    with open(path, 'r') as f:
        for line in f:
            if line.strip():
                values = line.strip().split('\t')
                if len(values) >= 3:
                    try:
                        contig = str(values[0])
                        position = int(values[1])
                        depth = float(values[2])
                        depth_list.append({
                            'contig': contig,
                            'position': position,
                            'depth': depth
                        })
                    except (ValueError, IndexError):
                        continue
    return depth_list

# Load depth and gene annotation
depth_list = load_depth(depth_path)

# Read marker gene annotation tsv (expects columns: gene_id, start, end)
def load_gene_file(path):
    genes = []
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line_num, line in enumerate(f, 2):
            if line.strip():
                values = line.strip().split('\t')
                if len(values) != len(header):
                    # Silent: skip malformed lines
                    pass
                    continue
                gene_dict = dict(zip(header, values))
                genes.append(gene_dict)
                if line_num <= 3:
                    pass
    # Silent
    return genes

gene_list = load_gene_file(gene_path)

gene_stats = []
phage_stats = []

# Marker gene statistics
for i, gr in enumerate(gene_list):
    if 'gene_id' not in gr:
        print("Error: gene_id not found in record {}: {}".format(i, gr))
        continue
    if 'start' not in gr:
        print("Error: start not found in record {}: {}".format(i, gr))
        continue
    if 'end' not in gr:
        print("Error: end not found in record {}: {}".format(i, gr))
        continue
        
    gid    = gr['gene_id']
    chrom  = '_'.join(gid.split('_')[:-1])
    start  = int(gr['start'])
    stop   = int(gr['end'])
    
    # Collect depths within the gene region
    seg_depths = []
    for depth_record in depth_list:
        if (depth_record['contig'] == str(chrom) and 
            depth_record['position'] >= start and 
            depth_record['position'] <= stop):
            seg_depths.append(depth_record['depth'])
    
    tot    = sum(seg_depths)
    length = stop - start + 1
    per    = tot/length if length > 0 else 0
    med    = statistics.median(seg_depths) if seg_depths else 0
    gene_stats.append({
        'gene_id':       gid,
        'Total_Counts':  tot,
        'Ave_Counts':    per,
        'Median_Depth':  med,
        'Region_Length': length
    })

# Phage region statistics (based on new phage_info.tsv format)
for pr in phage_list:
    contig = pr['contig']
    st = pr['start']
    sp = pr['end']

    # Collect depths within the phage region
    seg_depths = []
    for depth_record in depth_list:
        if (depth_record['contig'] == str(contig) and 
            depth_record['position'] >= st and 
            depth_record['position'] <= sp):
            seg_depths.append(depth_record['depth'])

    tot    = sum(seg_depths)
    length = sp - st + 1
    per    = tot/length if length > 0 else 0
    med    = statistics.median(seg_depths) if seg_depths else 0
    phage_stats.append({
        'Phage_Id':      pr['phage_name'],
        'Contig':        contig,
        'Start':         st,
        'Stop':          sp,
        'Total_Counts':  tot,
        'Ave_Counts':    per,
        'Median_Depth':  med,
        'Region_Length': length
    })

# Write outputs
gene_out  = os.path.join(output_dir, "marker_gene_counts.tsv")
phage_out = os.path.join(output_dir, "phage_counts.tsv")

# Write marker gene stats
with open(gene_out, 'w') as f:
    if gene_stats:
        # Header
        headers = list(gene_stats[0].keys())
        f.write('\t'.join(headers) + '\n')
        # Rows
        for stat in gene_stats:
            values = [str(stat[header]) for header in headers]
            f.write('\t'.join(values) + '\n')

# Write phage stats
with open(phage_out, 'w') as f:
    if phage_stats:
        # Header
        headers = list(phage_stats[0].keys())
        f.write('\t'.join(headers) + '\n')
        # Rows
        for stat in phage_stats:
            values = [str(stat[header]) for header in headers]
            f.write('\t'.join(values) + '\n')

# Compute host marker gene depth median
# Median of Ave_Counts across marker genes
ave_counts_list = [stat['Ave_Counts'] for stat in gene_stats]
median_total = statistics.median(ave_counts_list) if ave_counts_list else 0

# Host name is derived from the host genome file basename (without extension)
host_basename = os.path.basename(host_genome_path)
host_name = os.path.splitext(host_basename)[0]

# Write host counts with header: Host, Median_of_MG
host_out = os.path.join(output_dir, "host_counts.tsv")
with open(host_out, 'w') as f:
    f.write('Host\tMedian_of_MG\n')
    f.write('{}\t{}\n'.format(host_name, median_total))

# Silent

