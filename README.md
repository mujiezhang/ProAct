<div align="center">
  <img src="https://github.com/user-attachments/assets/ffe1315d-9e5c-48cd-a712-a3765102a4b5" alt="ProBord" width="450" />
</div>

# ProAct: **Pro**virus **Act**ivity Detector ✨
根据宿主测序原始数据评估其中原噬菌体的活跃度。

## 💠 目录

- [ProAct检测原理](#-ProAct检测原理)
- [ProAct工作流程](#-ProAct工作流程)
- [依赖环境](#-依赖环境)
- [安装与环境准备](#-安装与环境准备)
- [使用说明](#-使用说明)
  - [输入参数](#-输入参数)
  - [命令示例](#-命令示例)
  - [输出文件](#-输出文件)
- [测试数据](#-测试数据)
- [Citation](#-citation)
- [Contact](#-contact)

---

# Introduction

## 💠 Principle
Using whole-genome sequencing (WGS) data, ProAct exploits the principle that a provirus in lysogeny shares the same copy number as its host, resulting in a Provirus-to-Host coverage ratio (PtoH) of 1, whereas transition toward lysis drives self-replication and elevates PtoH above 1.

<img width="803" height="197" alt="schematic" src="https://github.com/user-attachments/assets/05ddcefd-5bdb-4298-8e27-0fe7ee55f065" />

## 💠 Workflow of ProBord
`ProAct` 需要输入宿主参考基因组及其原始测序数据、原噬菌体起始位点，通过（1）质控过滤后的原始读段和参考基因组进行比对，生成覆盖深度数据；（2）分别计算每个marker gene区域的平均coverage，取中位值代表宿主coverage；计算原噬菌体区域的平均coverage代表噬菌体coverage；（3）计算PtoH，得到该宿主内原噬菌体的活跃度（通过PtoH表征）。

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
- Install miniconda and add channels (If already installed, please skip)
```
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels bioconda
conda config --add channels conda-forge
```
- Install dependencies
``` 
conda create -n proact python3 bbmap bwa samtools hmmer prodigal 
conda activate probord
pip install pandas biopython pysam
```
- Download ProAct from github
```
git clone https://github.com/mujiezhang/ProAct.git
```

## 💠 How to run

- ▶️ Command line options
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
**Note**: `ERR4552622_R1.fastq` and `ERR4552622_R2.fastq` are the raw reads; `GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna` is the host genome; `phage_info.tsv` contains location information for `phage p22`, `phage Fels-1` and `phage Fels-2`. We get these information from paper `Turkington C J R, Abadi N N, Edwards R A, et al. hafeZ: active prophage identification through read mapping[J]. bioRxiv, 2021: 2021.07. 21.453177.`

- run an example
```
python ProAct/proact_pipeline.py -g test-data/GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna -1 test-data/ERR4552622_R1.fastq -2 test-data/ERR4552622_R2.fastq -p test-data/phage_info.tsv -t 40 -o test-result
```

## Output files
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
5. `PtoH_results.tsv`: prediction result

A detailed overview of `PtoH_results.tsv`:

| Host | Phage_Id | Contig | Start | Stop | Total_Counts | Ave_Counts | Median_of_MG | PtoH | Reads_depth_quality |
|------|----------|--------|-------|------|--------------|------------|--------------|------|---------------------|
| GCA_904129595.1_S.Tm_LT2p22_assembled_genomic | Sal_p22 | LR881463.1 | 1213987 | 1255756 | 149069351.0 | 3568.8137658606656 | 217.32838283828383 | 16.421296285613348 | high |
| GCA_904129595.1_S.Tm_LT2p22_assembled_genomic | Sal_Fels-1 | LR881463.1 | 1849458 | 1892188 | 39795679.0 | 931.3069902412768 | 217.32838283828383 | 4.285252474060286 | high |
| GCA_904129595.1_S.Tm_LT2p22_assembled_genomic | Sal_Fels-2 | LR881463.1 | 3731215 | 3764954 | 9870485.0 | 292.54549496147007 | 217.32838283828383 | 1.3460988902639377 | high |

- **`sample.depth`**：参考基因组上每个位点的测序深度（chromosome, position, depth）。
- **`sample_phage_info.txt`**：根据用户输入信息汇总（sample, genome, contig, start-end, phage_id）。
- **`MG/`**：GTDB-Tk 标注结果目录，包含：
  - `genome.tsv`：marker gene 注释表
  - `gtdbtk.json`：注释参数与版本信息
  - `gtdbtk.log` / `gtdbtk.warnings.log`：运行日志与警告
  - `identify/…`：HMM 比对与中间文件
- **`counts/`**：深度统计结果目录，包含：
  - `marker_gene_counts.tsv`：宿主基因组中每个 marker gene 的信息（Gene Id, Total_Counts, Per_Counts, Median_Depth, Region_Length）
  - `phage_counts.tsv`：每个噬菌体的信息（Phage_Id, Chromosome, Start, Stop, Total_Counts, Per_Counts, Median_Depth, Region_Length）
  - `host_counts.tsv`：宿主基因组中所有 marker gene 平均深度的中位值（Sample_ID, Median_of_MG）
- **`PtoH.tsv`**：最终输出，将 `phage_counts.tsv`的`Per_Counts` 除以 `host_counts.tsv`的`Median_of_MG` 得到 PtoH 值，并附加质量标签（`high`/`low`）。如果phage_Per_Counts > 10 且 host_Median_of_MG > 10，则为'high'；否则为'low'；如果PtoH ≥ 1.5，则活性判定为"active"，若PtoH < 1.5，则活性判定为"inactive".

  💡 **Note 2:** 
  如果phage_Per_Counts < 10 或 host_Median_of_MG < 10，则质量标签为'low'。表明测序深度不足，容易导致PtoH偏差较大。

# 💠 Citation
......

# 📬 Contact
```
# Mujie Zhang
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University
# Email: zhangmujie@sjtu.edu.cn
```
