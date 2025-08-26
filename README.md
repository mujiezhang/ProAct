<img width="396" height="366" alt="logo" src="https://github.com/user-attachments/assets/ffe1315d-9e5c-48cd-a712-a3765102a4b5" />

## ProAct:Provirus Activity Detector🦸
根据宿主测序原始数据评估其中原噬菌体的活跃度。

## 💠 目录

- [ProAct检测原理](#ProAct检测原理)
- [ProAct工作流程](#ProAct工作流程)
- [目录结构](#目录结构)
- [依赖环境](#依赖环境)
- [安装与环境准备](#安装与环境准备)
- [使用说明](#使用说明)
  - [输入参数](#输入参数)
  - [命令示例](#命令示例)
  - [输出文件](#输出文件)
- [测试数据](#测试数据)
- [引用](#引用)
- [贡献与反馈](#贡献与反馈)

---

## 💠 ProAct检测原理
Using whole-genome sequencing (WGS) data, ProAct exploits the principle that a provirus in lysogeny shares the same copy number as its host, resulting in a Provirus-to-Host coverage ratio (PtoH) of 1, whereas transition toward lysis drives self-replication and elevates PtoH above 1.

<img width="803" height="197" alt="schematic" src="https://github.com/user-attachments/assets/05ddcefd-5bdb-4298-8e27-0fe7ee55f065" />

## 💠 ProAct工作流程
`ProAct` 需要输入宿主参考基因组及其原始测序数据、原噬菌体起始位点，通过（1）质控过滤后的原始读段和参考基因组进行比对，生成覆盖深度数据；（2）分别计算每个marker gene区域的平均coverage，取中位值代表宿主coverage；计算原噬菌体区域的平均coverage代表噬菌体coverage；（3）计算PtoH，得到该宿主内原噬菌体的活跃度（通过PtoH表征）。

<img width="787" height="199" alt="workflow" src="https://github.com/user-attachments/assets/6c22bc29-d1eb-40cd-ad99-8127762a3adf" />

## 💠 目录结构

```plain
ProAct/
├── 1_form_depth.sh          # 质控过滤后的原始读段和参考基因组进行比对，生成 .depth 文件
├── 2_form_MG.sh             # GTDB-Tk 注释 marker genes
├── 3.get_MG_phage_counts.py # 计算每个 marker gene区域的平均coverage，取中位值代表宿主coverage；计算原噬菌体区域的平均coverage代表噬菌体coverage
├── 4.calculate_PtoH.py      # 计算 PtoH 指标
├── run_pipeline.sh         # 串联以上四步的脚本
├── test/                    # 测试数据目录
│   ├── GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna #宿主基因组
│   ├── ERR4552622_R1.fastq #宿主测序 R1 数据
│   └── ERR4552622_R2.fastq #宿主测序 R2 数据
└── README.md               # 本文档
```

## 💠 依赖环境

- **操作系统**: Linux / macOS
- **语言与工具**:
  - Bash >=4.0
  - Python 3.8+
  - [Conda](https://docs.conda.io/)
  - [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/)
  - BWA, Samtools
  - GTDB-Tk 2.x
  - Python 库: `pandas`, `biopython`, `pysam`

## 💠 安装与环境准备

你可以使用 `conda install` 命令逐个安装所有依赖：

```bash
# 创建并激活一个新环境
conda create -n ProAct bash python=3.8 -y
conda activate ProAct

# 安装脚本所需工具
conda install -c conda-forge -c bioconda bbmap bwa samtools gtdbtk -y

# 安装 Python 库
conda install pandas biopython pysam -y
```

## 💠 使用说明

### ▶️ 输入参数

| 参数   | 说明                         | 示例                            |
| ---- | -------------------------- | ----------------------------- |
| `-g` | 宿主基因组 FASTA 文件             | `genome.fna`                  |
| `-1` | 测序数据 R1 FASTQ              | `reads_R1.fastq`              |
| `-2` | 测序数据 R2 FASTQ              | `reads_R2.fastq`              |
| `-r` | 噬菌体区段：`contig,起始-终止`，可多次指定 | `-r CP024621.1,292609-333351` |
| `-o` | 输出目录                       | `results/`                    |
| `-t` | 比对线程数，默认 8                 | `-t 16`                       |
| `-c` | GTDB-Tk CPU 数，默认 4         | `-c 8`                        |
| `-e` | Conda 环境名，可选               | `-e ProAct`             |

💡 **Note 1:** 
  （1）测序数据需要为双端测序；
  （2）宿主基因组必须包含噬菌体所在的contig区域
  
### ▶️ 命令示例

```bash
bash run_pipeline.sh -g genome.fna -1 reads_R1.fastq -2 reads_R2.fastq -r contig1,1-5000 -r contig2,10000-20000 -o results/ -t 16 -c 8
```

### ▶️ 输出文件

运行完成后，`results/` 下会包含：

```
results/
├── sample.depth
├── sample_phage_info.txt
├── MG/
│   ├── genome.tsv
│   ├── gtdbtk.json
│   ├── gtdbtk.log
│   ├── gtdbtk.warnings.log
│   └──identify
│      └──...
├── counts/
│   ├── marker_gene_counts.tsv
│   ├── phage_counts.tsv
│   └── host_counts.tsv
└── PtoH.tsv
```
#### 🔹 输出结果说明
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

## 💠 测试数据
1.通过wget命令从zenodo下载测试数据（ERR4552622_R1.fastq, ERR4552622_R2.fastq,GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna）
```bash
# R1
wget -c 'https://zenodo.org/records/16909693/files/ERR4552622_R1.fastq?download=1'
# R2
wget -c 'https://zenodo.org/records/16909693/files/ERR4552622_R2.fastq?download=1'
# 参考基因组 .fna
wget -c 'https://zenodo.org/records/16909693/files/GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna?download=1'
```
2.激活ProAct环境，并运行run_pipeline.sh
```bash
conda activate ProAct
cd ProAct/
bash run_pipeline.sh -g test/GCA_904129595.1_S.Tm_LT2p22_assembled_genomic.fna -1 test/ERR4552622_R1.fastq -2 test/ERR4552622_R2.fastq -r LR881463.1,1213987-1255756 -r LR881463.1,1849458-1892188 -r LR881463.1,3731215-3764954 -o result/ -t 10 -c 8
```
3.结果文件：
3.1 result/目录下文件如下所示：

<img width="671" height="467" alt="result" src="https://github.com/user-attachments/assets/a4d012f2-4c08-4560-8a58-41b573159fac" />

3.2 PtoH.tsv文件为判定结果文件，内容如图所示：



## 💠 引用
...

## 💠 贡献与反馈

欢迎提出 issue 或 pull request，共同完善项目。
