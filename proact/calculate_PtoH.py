#!/usr/bin/env python3
"""
Calculate PtoH (Provirus-to-Host coverage ratio).

For each phage record: PtoH = Ave_Counts / Median_of_MG
Reads_depth_quality = "low" if (Ave_Counts < 10 or Median_of_MG < 10) else "high"
Predicted_activity = "active" if PtoH >= 1.5 else "inactive"
"""

import sys
import os
import argparse
import pandas as pd


def compute_ptoh(phage_file, host_file, output_file):
    """Compute PtoH from phage and host depth statistics."""
    phage_df = pd.read_csv(phage_file, sep="\t")
    host_df = pd.read_csv(host_file, sep="\t")

    host_name = host_df.loc[0, "Host"]
    median_mg = host_df.loc[0, "Median_of_MG"]

    phage_df["Host"] = host_name
    phage_df["Median_of_MG"] = median_mg
    phage_df["PtoH"] = phage_df["Ave_Counts"] / phage_df["Median_of_MG"]

    phage_df["Reads_depth_quality"] = phage_df.apply(
        lambda r: "low" if (r["Ave_Counts"] < 10 or r["Median_of_MG"] < 10) else "high",
        axis=1
    )

    phage_df["Predicted_activity"] = phage_df["PtoH"].apply(
        lambda x: "active" if x >= 1.5 else "inactive"
    )

    cols = [
        "Host", "Phage_Id", "Contig", "Start", "Stop",
        "Total_Counts", "Ave_Counts", "Median_of_MG", "PtoH",
        "Predicted_activity", "Reads_depth_quality",
    ]
    out_df = phage_df[cols]

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    out_df.to_csv(output_file, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        prog='proact-calculate-ptoh',
        description='Calculate PtoH from phage and host depth statistics',
    )
    parser.add_argument('phage_file', help='phage_counts.tsv')
    parser.add_argument('host_file', help='host_counts.tsv')
    parser.add_argument('output_file', help='Output TSV file')

    args = parser.parse_args()
    compute_ptoh(args.phage_file, args.host_file, args.output_file)


if __name__ == "__main__":
    main()
