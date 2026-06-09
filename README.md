# PHB_gtdb — GTDB 数据库中 PHB 塑料降解基因 (PhaZ) 的系统发育分析

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 简介

本项目从 **GTDB** (Genome Taxonomy Database) Release R232 中系统鉴定具有 **PHB (聚羟基丁酸酯)** 塑料降解能力的细菌，并进行完整的生物信息学分析。

### 核心目标

- **基因搜索**: 使用 **14 条 NCBI 验证的 PhaZ 参考序列** (Pyrodigal + DIAMOND)，搜索 165,468 个细菌代表基因组
- **古菌验证**: 同步搜索古菌基因组 (~34,455 个)，验证"古菌不携带 phaZ 基因"的假说
- **CD-HIT 去冗余**: c95 去真冗余；c70 仅用于概览热图，不作为系统发育证据
- **结构域分型**: 按 Pfam 结构域 + 参考序列类型分为 5 个亚型 (胞外2型 + 胞内2型 + 新型)
- **同型独立建树**: 各亚型独立 MAFFT + trimAl + IQ-TREE，不混合不同类别
- **可视化**: 系统发育树 + c70 概览热图 + 门分布

### PhaZ 结构域分型

PhaZ 是 identity 跨度 30%–99% 的 α/β 水解酶超家族。按 Pfam 结构域分为:

| 亚型 | 序列数 | 信号肽 | Lipase box | 代表物种 | 定位 |
|------|:---:|:---:|:---:|---|---|
| **extracellular** | 596 | ✓ | ✓ (G-X-S-X-G) | *P. stutzeri*, *C. testosteroni* | 胞外分泌 |
| **extracellular_lemoignei** | 434 | ✓ | ✓ | *P. lemoignei* (7个旁系同源) | 胞外分泌 |
| **intracellular** | 4,693 | ✗ | ✗ | *C. necator* H16 PhaZ1/2/5 | 胞内 |
| **ralstonia** | 2,948 | ✗ | ✗ | *R. pickettii* | 胞内变体 |
| **bacillus_type** | 60 | ✗ | ✓ (G-W-S-M-G) | *Bacillus* sp. CDB3 | 胞内 (革兰氏+) |

> 详细分类依据见 [`docs/PIPELINE.md`](docs/PIPELINE.md)

### CD-HIT 去冗余实验

| 阈值 | 保留率 | 评价 |
|:---:|:---:|---|
| **c95** | **99.6%** | ✓ 仅去真冗余, 保留全部功能多样性 |
| c90 | 98.4% | 几乎无变化 |
| c70 | 73.9% | 开始丢失亚家族多样性 |
| c60 | 50.2% | ✗ 不同亚家族被强制合并 |
| c50 | ~28% | ✗ 功能信号严重失真 |

> **原则**: c95 去真冗余 + 结构域分型 + 同型各自建树。c70 仅用于概览热图，**绝不作为系统发育证据**。

### 参考序列

> ⚠️ **重要**: 所有 PhaZ 参考序列均通过 NCBI E-utilities API 实时查询验证，确保每一条 accession 对应正确的 PHB/PHA depolymerase 蛋白。详见 [`docs/PHAZ_REFERENCES.md`](docs/PHAZ_REFERENCES.md)。

**14 条参考序列** 覆盖 3 个门、8 个属，包含胞外型 (extracellular) 和胞内型 (intracellular)：

| 类型 | 代表物种 | 序列数 | 代表性文献 |
|------|---------|--------|-----------|
| 胞外型 (extracellular) | *Stutzerimonas stutzeri*, *Comamonas testosteroni*, *Paucimonas lemoignei* 等 | 8 | Ohura et al., 1999; Jendrossek et al., 1995 |
| 胞内型 (intracellular) | *Cupriavidus necator* H16 | 3 | Saegusa et al., 2001 |
| 新型 (new type) | *Bacillus* sp. CDB3, *Ralstonia pickettii* | 3 | Tseng et al., 2006 |

### 多拷贝基因策略

部分基因组含有多个 PhaZ 拷贝 (最高 12 个)。当前策略: **保留所有拷贝**, 在系统发育树中标注基因组 ID。拷贝数分布统计将作为功能注释的一部分。

### IQ-TREE 建树参数

| 亚型 | 序列 | 模型 | Bootstrap | 策略 |
|---|---|---|---|---|
| bacillus_type | 57 | MFP | B1000 | ✅ 完成 |
| extracellular_lemoignei | 412 | MFP | B1000 | ✅ 完成 |
| extracellular | 509 | MFP | B1000 | ✅ 完成 |
| ralstonia | 2,776 | LG+F+G4 | — | ML only |
| intracellular | 4,424 | LG+F+G4 | — | ML only |

> 大组跳过 ModelFinder（400+ 模型测试 prohibitive），使用通用 LG+F+G4 替代矩阵 + `-T AUTO` 线程。

### 关键结果 (截至 2026-06-09)

| 指标 | 值 |
|---|---|
| 搜索基因组 | 165,468 细菌 + 10,122 古菌 |
| DIAMOND 命中 | 7,068 细菌 + 2 古菌 |
| HMMer 验证率 | 80.7% (6,033/7,478) |
| **最终验证集** | **6,532 条 PhaZ** (15 细菌门, **0 古菌**) |
| 亚型 | intracellular (56%), ralstonia (34%), extracellular (6%), lemoignei (4%), bacillus (<1%) |
| Lipase box | 胞外型 28-53% vs 胞内型 ~10% → 分型验证通过 |
| 古菌结论 | 2/2 假阳性 (34-44 aa 片段) → **古菌不携带 phaZ** |

详细统计: [`RESULTS.md`](docs/RESULTS.md) | 方法: [`METHODS.md`](docs/METHODS.md)

### 改进计划

已完成: HMMer 验证 ✅ | Lipase box 检测 ✅ | 古菌验证 ✅ | 长度过滤 ✅ | 拷贝数统计 ✅
待完成: SignalP | InterProScan | 可视化热图 | 联合树

详见 [`docs/REVIEW.md`](docs/REVIEW.md) — 21 项改进路线图

---

## 项目结构

```
PHB_gtdb/
├── README.md
├── docs/
│   ├── PIPELINE.md                   # 完整分析管线文档
│   ├── METHODS.md                    # 方法与结果详细报告 (论文用)
│   ├── PHAZ_REFERENCES.md            # PhaZ 参考序列 (14条, NCBI验证)
│   ├── RESULTS.md                    # 最终统计结果 (门矩阵, 亚型, Top属)
│   ├── SEARCH_REPORT.md              # 搜索结果报告
│   └── REVIEW.md                     # 管线评审与改进计划
├── scripts/
│   ├── config.py                     # 全局配置
│   ├── utils.py                      # 工具函数
│   ├── 01_phb_search.py              # Step 1: 细菌 PhaZ 搜索
│   ├── 01b_archaea_search.py         # Step 1b: 古菌 PhaZ 搜索
│   ├── 02_extract_sequences.py       # Step 2: 提取 PhaZ 蛋白序列
│   ├── 03_msa.py                     # Step 3: 分型 + MAFFT + trimAl
│   ├── 04_phylogeny.py              # Step 4: IQ-TREE 建树
│   ├── 05_annotation.py             # Step 5: 功能注释
│   └── 06_visualization.py          # Step 6: 可视化
├── data/processed/                   # 中间结果
├── results/                          # 输出 (树文件等)
├── monitor_search.sh                 # 实时搜索监控
├── monitor.sh                        # 通用分析监控
├── download_monitor.sh               # 下载进度监控
└── extract_monitor.sh                # 解压进度监控
```

## 快速开始

### 环境要求

- **Python** ≥ 3.12
- **Conda** (推荐) 或 pip + 系统工具
- **生物信息学工具**: DIAMOND, MAFFT, trimAl, IQ-TREE

### 安装

```bash
# 创建 conda 环境
conda create -n phb_gtdb python=3.12
conda activate phb_gtdb

# 安装 Python 依赖
pip install -r requirements.txt

# 安装生物信息学工具
conda install -c bioconda -c conda-forge \
    diamond mafft trimal iqtree
```

### 运行

```bash
# === 细菌搜索 ===
# 测试模式 (100 个基因组)
python scripts/01_phb_search.py --threads 30 --limit 100

# 全量搜索 (165,468 个基因组, ~38h)
nohup conda run -n phb_gtdb python scripts/01_phb_search.py --threads 30 &

# === 古菌搜索 ===
# 全量搜索 (~34,455 个基因组)
nohup conda run -n phb_gtdb python scripts/01b_archaea_search.py --threads 30 &

# === 实时监控 ===
bash monitor_search.sh                # 一次性状态
bash monitor_search.sh --watch 15     # 每15秒自动刷新
bash monitor_search.sh --mode archaea # 指定搜索模式

# === 后续步骤 ===
python scripts/02_extract_sequences.py
python scripts/03_msa.py --threads 30
python scripts/04_phylogeny.py --threads 30
python scripts/05_annotation.py --threads 30
python scripts/06_visualization.py
```

## 搜索方法

### 基因预测: Pyrodigal

使用 [Pyrodigal](https://github.com/althonos/pyrodigal) (Prodigal 的 Python 绑定) 进行**进程内**基因预测，避免 subprocess 开销：

```python
orf_finder = pyrodigal.GeneFinder(meta=True)
genes = orf_finder.find_genes(dna_seq)
for gene in genes:
    protein_seq = gene.sequence()
```

### 同源搜索: DIAMOND blastp

DIAMOND 比传统 BLASTP 快 100-1000x：

```
diamond blastp -q proteins.faa -d phaz_db.dmnd -o output.tsv \
    --outfmt 6 qseqid sseqid pident evalue \
    -e 1e-10 --id 30 --query-cover 50 --max-target-seqs 5
```

### 筛选标准

| 参数 | 值 | 说明 |
|------|-----|------|
| E-value | ≤ 1e-10 | 严格的统计显著性 |
| Identity | ≥ 30% | 允许种间变异 |
| Query coverage | ≥ 50% | 确保覆盖催化域 |

### 性能优化

- **chunksize=1**: 最优负载均衡（慢基因组不会阻塞批次）
- **增量保存**: 每 500 个基因组保存一次，防止数据丢失
- **断点续传**: 重启后自动跳过已处理基因组
- **30 线程并行**: 充分利用多核 CPU

## 数据源

GTDB 代表基因组数据库 (Release R232) 位于服务器:
```
/home/data/haoyu/GTDB/
├── gtdb_genomes_reps_r232/database/   # 199,923 个基因组 (GCA+GCF)
├── metadata/bac120_metadata_r232.tsv  # 质量元数据
├── taxonomy/bac120_taxonomy_r232.tsv  # 分类学信息
└── GTDB_tree/bac120_r232.tree         # 物种系统发育树
```

## 引用

- GTDB: Parks et al. (2022) *Nature Biotechnology* 40:1273-1281
- Pyrodigal: Larralde (2022) *JOSS* 7(72):4296
- DIAMOND: Buchfink et al. (2021) *Nature Methods* 18:366-368
- MAFFT: Katoh & Standley (2013) *MBE* 30:772-780
- IQ-TREE: Minh et al. (2020) *MBE* 37:1530-1534
- PhaZ 参考序列: 详见 [`docs/PHAZ_REFERENCES.md`](docs/PHAZ_REFERENCES.md)

## 许可证

MIT License
