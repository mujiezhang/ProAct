#!/usr/bin/env python3
"""
Script: calculate_PtoH.py

Inputs:
  1) phage_counts.tsv (per-phage region depth)
  2) host_counts.tsv (median depth of marker genes)

Compute:
  For each phage record: PtoH = Ave_Counts / Median_of_MG
  Reads_depth_quality = "low" if (Ave_Counts < 10 or Median_of_MG < 10) else "high"

Output columns:
  Host, Phage_Id, Contig, Start, Stop, Total_Counts, Ave_Counts, Median_of_MG, PtoH, Reads_depth_quality
"""

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <phage_counts.tsv> <host_counts.tsv> <output.tsv>")
        sys.exit(1)

    phage_file, host_file, output_file = sys.argv[1:]

    # Read inputs
    phage_df = pd.read_csv(phage_file, sep="\t")
    host_df = pd.read_csv(host_file, sep="\t")

    # Use the first host record
    host_name   = host_df.loc[0, "Host"]
    median_mg   = host_df.loc[0, "Median_of_MG"]

    # Broadcast to each phage record
    phage_df["Host"]          = host_name
    phage_df["Median_of_MG"]  = median_mg
    phage_df["PtoH"]          = phage_df["Ave_Counts"] / phage_df["Median_of_MG"]

    # Reads depth quality
    phage_df["Reads_depth_quality"] = phage_df.apply(
        lambda r: "low" if (r["Ave_Counts"] < 10 or r["Median_of_MG"] < 10) else "high",
        axis=1
    )

    # Predicted activity
    phage_df["Predicted_activity"] = phage_df["PtoH"].apply(
        lambda x: "active" if x >= 1.5 else "inactive"
    )

    # Reorder columns
    cols = [
        "Host", "Phage_Id", "Contig", "Start", "Stop",
        "Total_Counts", "Ave_Counts", "Median_of_MG", "PtoH", "Predicted_activity", "Reads_depth_quality"
    ]
    out_df = phage_df[cols]

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    out_df.to_csv(output_file, sep="\t", index=False)
    # Silent

if __name__ == "__main__":
    main()

