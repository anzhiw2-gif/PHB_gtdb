# PHB_gtdb — GTDB 数据库中 PHB 塑料降解基因的系统发育与功能分析

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 简介

本项目从 **GTDB** (Genome Taxonomy Database) 公共数据库中系统鉴定具有 **PHB (聚羟基丁酸酯)** 塑料降解能力的古菌和细菌，并进行完整的生物信息学分析。

### 分析内容
- **基因搜索**: 使用 HMM 模型在 143,614 个 GTDB 代表基因组中搜索 PHB 降解/合成相关基因
- **序列比对**: MAFFT 多序列比对 + trimAl 修剪
- **系统发育分析**: IQ-TREE 最大似然树 (基因树 + 物种树)
- **功能注释**: eggNOG-mapper (KEGG / COG / GO)
- **可视化**: 论文级别图表 (门分布、热图、通路完整性等)

### 分析基因
| 基因 | 功能 | Pfam ID |
|------|------|---------|
| **phaZ** | PHB depolymerase (降解) | PF14556 |
| **phaC** | PHA synthase Class I (合成) | PF07167 |
| **phaC2** | PHA synthase Class II | PF06863 |
| **phaA** | β-ketothiolase | K00626 |
| **phaB** | Acetoacetyl-CoA reductase | K00023 |
| **phaP** | Phasin (PHA granule protein) | PF09650 |

## 项目结构

```
PHB_gtdb/
├── main_pipeline.sh              # 一键运行完整流程
├── README.md
├── requirements.txt
├── .gitignore
├── scripts/
│   ├── config.py                 # 全局配置
│   ├── utils.py                  # 工具函数
│   ├── 01_phb_search.py          # Step 1: HMM 搜索
│   ├── 02_extract_sequences.py   # Step 2: 序列提取
│   ├── 03_msa.py                 # Step 3: 多序列比对
│   ├── 04_phylogeny.py           # Step 4: 系统发育分析
│   ├── 05_annotation.py          # Step 5: 功能注释
│   └── 06_visualization.py       # Step 6: 可视化
├── data/
│   ├── external/phb_hmm/         # PHB 相关 Pfam HMM 模型
│   └── processed/                # 中间结果
├── results/
│   ├── figures/                  # 论文图表 (PDF)
│   ├── tables/                   # 统计表格 (TSV)
│   └── logs/                     # 运行日志
└── notebooks/                    # Jupyter Notebooks
```

## 快速开始

### 环境要求

- **Python** ≥ 3.10
- **系统工具**: MAFFT, trimAl, FastTree, BLAST+, HMMER, Prodigal
- **推荐**: IQ-TREE (更精确的 ML 树)

### 安装

```bash
# 安装 Python 依赖
pip install -r requirements.txt

# 安装生物信息学工具 (如果未安装)
conda install -c bioconda -c conda-forge \
    mafft trimal fasttree iqtree hmmer prodigal blast eggnog-mapper
```

### 运行

```bash
# 一键运行完整流程 (30 线程)
bash main_pipeline.sh --threads 30

# 或分步运行
python scripts/01_phb_search.py --threads 30     # HMM 搜索
python scripts/02_extract_sequences.py            # 序列提取
python scripts/03_msa.py --threads 30             # 多序列比对
python scripts/04_phylogeny.py --threads 30       # 系统发育
python scripts/05_annotation.py --threads 30      # 功能注释
python scripts/06_visualization.py                # 可视化

# 测试模式 (仅处理前 100 个基因组)
python scripts/01_phb_search.py --threads 30 --limit 100
```

### 数据源

GTDB 代表基因组数据库 (Release R226) 位于服务器:
```
/home/data/haoyu/GTDB/
├── gtdb_genomes_reps_r226/database/    # 143,614 个基因组
├── metadata/bac120_metadata_r226.tsv   # 质量元数据
├── taxonomy/bac120_taxonomy_r226.tsv   # 分类学信息
└── GTDB_tree/bac120_r226.tree          # 物种系统发育树
```

## 结果示例

![PHB Summary Dashboard](results/figures/phb_summary_dashboard.pdf)

## 引用

- GTDB: Parks et al. (2022) *Nature Biotechnology*
- MAFFT: Katoh & Standley (2013) *Molecular Biology and Evolution*
- IQ-TREE: Nguyen et al. (2015) *Molecular Biology and Evolution*
- eggNOG-mapper: Cantalapiedra et al. (2021) *Molecular Biology and Evolution*

## 许可证

MIT License
