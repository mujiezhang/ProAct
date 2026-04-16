#!/usr/bin/env python3
"""
Compute depth statistics for marker genes and phage regions.

Provides both a functional API (compute_depth_stats) and a CLI entry point.
"""
import os
import sys
import statistics
import argparse


def load_phage_info(path):
    """Load phage info from a tab-delimited file (name, contig, start, end)."""
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
                'end': end,
            })
    return phage_list


def load_depth(path):
    """Load depth data from a tab-delimited file (contig, position, depth)."""
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
                            'depth': depth,
                        })
                    except (ValueError, IndexError):
                        continue
    return depth_list


def load_gene_file(path):
    """Read marker gene annotation TSV (expects columns: gene_id, start, end)."""
    genes = []
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            if line.strip():
                values = line.strip().split('\t')
                if len(values) != len(header):
                    continue
                gene_dict = dict(zip(header, values))
                genes.append(gene_dict)
    return genes


def compute_depth_stats(phage_info_path, host_genome_path, depth_path,
                        gene_path, output_dir):
    """Compute depth statistics for marker genes and phage regions."""
    phage_list = load_phage_info(phage_info_path)
    depth_list = load_depth(depth_path)
    gene_list = load_gene_file(gene_path)

    os.makedirs(output_dir, exist_ok=True)

    gene_stats = []
    for i, gr in enumerate(gene_list):
        if 'gene_id' not in gr or 'start' not in gr or 'end' not in gr:
            print("Error: missing field in record {}: {}".format(i, gr))
            continue

        gid = gr['gene_id']
        chrom = '_'.join(gid.split('_')[:-1])
        start = int(gr['start'])
        stop = int(gr['end'])

        seg_depths = [
            d['depth'] for d in depth_list
            if d['contig'] == str(chrom)
            and d['position'] >= start
            and d['position'] <= stop
        ]

        tot = sum(seg_depths)
        length = stop - start + 1
        per = tot / length if length > 0 else 0
        med = statistics.median(seg_depths) if seg_depths else 0
        gene_stats.append({
            'gene_id': gid,
            'Total_Counts': tot,
            'Ave_Counts': per,
            'Median_Depth': med,
            'Region_Length': length,
        })

    phage_stats = []
    for pr in phage_list:
        contig = pr['contig']
        st = pr['start']
        sp = pr['end']

        seg_depths = [
            d['depth'] for d in depth_list
            if d['contig'] == str(contig)
            and d['position'] >= st
            and d['position'] <= sp
        ]

        tot = sum(seg_depths)
        length = sp - st + 1
        per = tot / length if length > 0 else 0
        med = statistics.median(seg_depths) if seg_depths else 0
        phage_stats.append({
            'Phage_Id': pr['phage_name'],
            'Contig': contig,
            'Start': st,
            'Stop': sp,
            'Total_Counts': tot,
            'Ave_Counts': per,
            'Median_Depth': med,
            'Region_Length': length,
        })

    # Write marker gene stats
    gene_out = os.path.join(output_dir, "marker_gene_counts.tsv")
    with open(gene_out, 'w') as f:
        if gene_stats:
            headers = list(gene_stats[0].keys())
            f.write('\t'.join(headers) + '\n')
            for stat in gene_stats:
                values = [str(stat[h]) for h in headers]
                f.write('\t'.join(values) + '\n')

    # Write phage stats
    phage_out = os.path.join(output_dir, "phage_counts.tsv")
    with open(phage_out, 'w') as f:
        if phage_stats:
            headers = list(phage_stats[0].keys())
            f.write('\t'.join(headers) + '\n')
            for stat in phage_stats:
                values = [str(stat[h]) for h in headers]
                f.write('\t'.join(values) + '\n')

    # Host marker gene depth median
    ave_counts_list = [stat['Ave_Counts'] for stat in gene_stats]
    median_total = statistics.median(ave_counts_list) if ave_counts_list else 0

    host_basename = os.path.basename(host_genome_path)
    host_name = os.path.splitext(host_basename)[0]

    host_out = os.path.join(output_dir, "host_counts.tsv")
    with open(host_out, 'w') as f:
        f.write('Host\tMedian_of_MG\n')
        f.write('{}\t{}\n'.format(host_name, median_total))


def main():
    parser = argparse.ArgumentParser(
        prog='proact-get-depth',
        description='Compute depth statistics for marker genes and phage regions',
    )
    parser.add_argument('phage_info', help='Phage info TSV file')
    parser.add_argument('host_genome', help='Host genome FASTA file')
    parser.add_argument('depth_file', help='Depth TSV file')
    parser.add_argument('marker_genes', help='Marker genes TSV file')
    parser.add_argument('output_dir', help='Output directory')

    args = parser.parse_args()
    compute_depth_stats(
        phage_info_path=args.phage_info,
        host_genome_path=args.host_genome,
        depth_path=args.depth_file,
        gene_path=args.marker_genes,
        output_dir=args.output_dir,
    )


if __name__ == '__main__':
    main()
