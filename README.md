<div align="center">
  <img src="https://github.com/user-attachments/assets/ffe1315d-9e5c-48cd-a712-a3765102a4b5" alt="ProBord" width="450" />
</div>

# ProAct: **Pro**virus **Act**ivity Detector ✨
Assessing the activity of provirus based on the host's raw sequencing data.

## Table of contents
<!-- TOC -->
- [ProAct: Provirus Activity Detector](#proact-provirus-activity-detector-)
- [Introduction](#introduction)
  - [Principle](#-principle)
  - [Workflow of ProAct](#-workflow-of-proact)
- [Instructions](#instructions)
  - [Dependencies](#-dependencies)
  - [**Installation**](#-installation)
  - [**How to run**](#-how-to-run)
  - [**Output files**](#-output-files)
- [Citation](#-citation)
- [Contact](#-contact)

<!-- /TOC -->


---

# Introduction

## 💠 Principle
Using whole-genome sequencing (WGS) data, ProAct exploits the principle that a provirus in lysogeny shares the same copy number as its host, resulting in a Provirus-to-Host coverage ratio (PtoH) of 1, whereas transition toward lysis drives self-replication and elevates PtoH above 1.

<img width="803" height="197" alt="schematic" src="https://github.com/user-attachments/assets/05ddcefd-5bdb-4298-8e27-0fe7ee55f065" />

## 💠 Workflow of ProAct
ProAct requires the input of the host reference genome, its original sequencing data, and the start/end site of the prophage. It proceeds by (1) aligning the quality-controlled and filtered raw reads to the reference genome to generate coverage depth data; (2) calculating the average coverage for each marker gene region, taking the median value to represent the host coverage, and calculating the average coverage of the prophage region to represent the phage coverage; (3) computing PtoH to obtain the activity level of the prophage within the host (represented by PtoH).

Note: We used the identify module of GTDB-Tk to identify marker genes. To avoid downloading the entire GTDB-Tk database when packaging the ProAct workflow, we employed the `extract_gtdb_mg.py` script from `https://github.com/4less/extract_gtdb_mg` with minor modifications to adapt it to ProAct's requirements. This script was adapted from GTDB-Tk, and its identification results are consistent with those of GTDB-Tk.

<img width="787" height="199" alt="workflow" src="https://github.com/user-attachments/assets/6c22bc29-d1eb-40cd-ad99-8127762a3adf" />

# Instructions

## 💠 Dependencies
```
python3
bbmap
bwa
samtools
hmmer
prodigal
pandas
biopython
pysam
```

## 💠 Installation
- Install miniconda and add channels (**If already installed, please skip**)
```
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels bioconda
conda config --add channels conda-forge
```
- Install dependencies
``` 
conda create -n proact python bbmap bwa samtools hmmer prodigal 
conda activate probord
pip install pandas biopython pysam
```
- Download ProAct from github
```
git clone https://github.com/mujiezhang/ProAct.git
```

## 💠 How to run

- ▶️ Command line options - `python ProAct/proact_pipeline.py -h`:
```
usage: proact_pipeline.py [-h] -g GENOME -1 READ1 -2 READ2 -p PHAGE_INFO [-o OUTPUT_DIR] [-t THREADS]
                          [-m] [--keep-tmp]

ProAct: Provirus Activity Detector

options:
  -h, --help            show this help message and exit
  -g, --genome GENOME   Host genome FASTA file
  -1, --read1 READ1     Read1 FASTQ file
  -2, --read2 READ2     Read2 FASTQ file
  -p, --phage-info PHAGE_INFO
                        A tab-delimited file with columns: virla_name, host_contig, start, end
  -o, --output-dir OUTPUT_DIR
                        Output directory (default: proact_results)
  -t, --threads THREADS
                        Number of threads (default: 40)
  -m, --meta            Run prodigal in meta mode (for metagenomes)
  --keep-tmp            Keep temporary files

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
```

Download test data:

```
wget -c https://zenodo.org/records/16948956/files/test-data.tar.gz
tar -xzvf test-data.tar.gz
```

Files in test-data:

```
test-data
├── ERR4552622_R1.fastq 
├── ERR4552622_R2.fastq
├── GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna
└── phage_info.tsv
```
**Note**: The files `ERR4552622_R1.fastq` and `ERR4552622_R2.fastq` contain the raw sequencing reads; `GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna` corresponds to the host genome; and `phage_info.tsv` provides the genomic location information for phage p22, phage Fels-1, and phage Fels-2. These data were obtained from the study by Turkington et al. (_hafeZ: active prophage identification through read mapping. bioRxiv, 2021: 2021.07.21.453177_).

- run an example
```
python ProAct/proact_pipeline.py -g test-data/GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna -1 test-data/ERR4552622_R1.fastq -2 test-data/ERR4552622_R2.fastq -p test-data/phage_info.tsv -t 40 -o test-result
```

## 💠 Output files
In this example, the results of ProBord's analysis will be written to the `test-result` directory, which will look like this:
```
test-result/
├── host_counts.tsv
├── marker_gene_counts.tsv
├── marker_genes.tsv
├── phage_counts.tsv
└── PtoH_results.tsv
```
1. `host_counts.tsv`: median depth of host marker gene
2. `marker_gene_counts.tsv`: depth information of host marker gene
3. `marker_genes.tsv`: identified host marker gene
4. `phage_counts.tsv`: depth information of host marker gene
5. `PtoH_results.tsv`: **prediction result**

A detailed overview of `PtoH_results.tsv`:

| Host | Phage_Id | Contig | Start | Stop | Total_Counts | Ave_Counts | Median_of_MG | PtoH | Predicted_activity | Reads_depth_quality |
|------|----------|--------|-------|------|--------------|------------|--------------|------|---------------------|---------------------|
| GCA_904129595.1_S.Tm_LT2p22_assembled_genomic | Sal_p22 | LR881463.1 | 1213987 | 1255756 | 149069575.0 | 3568.8191285611683 | 217.32838283828383 | 16.421320961177727 | active | high |
| GCA_904129595.1_S.Tm_LT2p22_assembled_genomic | Sal_Fels-1 | LR881463.1 | 1849458 | 1892188 | 39795641.0 | 931.3061009571504 | 217.32838283828383 | 4.285248382168952 | active | high |
| GCA_904129595.1_S.Tm_LT2p22_assembled_genomic | Sal_Fels-2 | LR881463.1 | 3731215 | 3764954 | 9870485.0 | 292.54549496147007 | 217.32838283828383 | 1.3460988902639377 | inactive | high |

- **`Predicted_activity`**： If PtoH ≥ 1.5, the activity is determined as "active"; if PtoH < 1.5, the activity is determined as "inactive".
- `Reads_depth_quality`: If Ave_Counts > 10 and Median_of_MG > 10, it is classified as 'high'; otherwise, it is classified as 'low' (indicating insufficient sequencing depth, which may lead to significant bias in PtoH).


# 💠 Citation
......

# 📬 Contact
```
# Mujie Zhang
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University
# Email: zhangmujie@sjtu.edu.cn
```
