# PHB_gtdb: GTDB R232 中 PhaZ/PHB 降解基因的系统分析

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)
[![R/gRodon2](https://img.shields.io/badge/R-gRodon2-green.svg)](https://github.com/jlw-ecoevo/gRodon2)

本项目从 GTDB Release R232 代表基因组中系统鉴定 PHB/PHA depolymerase, 即 `phaZ` 编码的 PHB 降解相关蛋白，并进一步分析其分类分布、亚型组成、保守 lipase box、系统发育关系，以及同属内 `phaZ+` 与 `phaZ-` 菌株的预测最大生长速率差异。

项目运行主体在 T141 服务器：

```text
Server: T141, 10.16.1.141
Project path: /home/data/haoyu/PHB_gtdb
GTDB root: /home/data/haoyu/GTDB
GTDB release: R232
```

本仓库保留可复现脚本、方法文档、汇报材料、Nature 风格论文图和作图源数据。GTDB 原始基因组、完整中间结果和大体积外部数据库不随 GitHub 上传。

## 主要结论

| 结论 | 当前结果 |
|---|---:|
| 搜索背景 | 165,468 个细菌代表基因组，10,122 个古菌基因组 |
| 细菌初始命中 | 7,068 个基因组，16,486 个 PhaZ 初始命中 |
| 最终高置信度 PhaZ | 6,532 条蛋白序列 |
| 覆盖范围 | 15 个细菌门，1,135 个属 |
| 古菌结果 | 初始 2 个短片段命中，最终 0 条高置信度 PhaZ |
| 主导类群 | Pseudomonadota 占主导 |
| 主导亚型 | 胞内型 PhaZ 为主 |
| 生长速率扩展分析 | gRodon2 同属 `phaZ+` vs `phaZ-`: 8,788 基因组, 899 属, **无显著差异** (p=0.459) |

核心解释：PhaZ 不是单一高度保守蛋白，而是 alpha/beta hydrolase 相关的功能多样酶家族。因此本项目采用“先高置信度筛选，再按亚型独立分析”的路线，避免把胞外、胞内和新型 PhaZ 混在一棵树中造成错误解释。

## 分析路线

```text
GTDB R232 representative genomes
  -> Pyrodigal 预测 CDS 和蛋白
  -> DIAMOND 搜索 14 条验证 PhaZ 参考序列
  -> 提取候选 PhaZ
  -> CD-HIT c95 去近完全重复
  -> 长度过滤，去除短片段
  -> HMMER profile 验证 + 高置信度 DIAMOND 补充
  -> 根据参考序列、结构特征和 lipase box 划分亚型
  -> 各亚型独立 MAFFT/trimAl 比对
  -> IQ-TREE 或 FastTree 建树
  -> 分类学分布、亚型组成和图形汇总
```

扩展分析：

```text
GTDB genomes in phaZ-containing genera
  -> Pyrodigal 预测 CDS 和蛋白
  -> Pfam ribosomal protein HMM + HMMER 识别核糖体蛋白
  -> 核糖体蛋白作为 highly expressed genes
  -> gRodon2 predictGrowth 预测最大生长速率
  -> 同属内比较 phaZ+ 与 phaZ-
```

## PhaZ 亚型命名

专家建议后，推荐在后续报告中使用“定位大类 + 谱系后缀”的命名方式，保留生物学解释同时避免把胞外/胞内误解为唯一分类依据。

| 推荐名称 | 原始名称 | 解释 |
|---|---|---|
| `extracellular.general` | `extracellular` | 一般胞外分泌型 PhaZ |
| `extracellular.lemoignei_like` | `extracellular_lemoignei` | *Paucimonas lemoignei* 相关胞外谱系 |
| `intracellular.cupriavidus_like` | `intracellular` | *Cupriavidus necator* H16 PhaZ1/2/5 相关胞内型 |
| `intracellular.ralstonia_like` | `ralstonia` | *Ralstonia* 相关胞内变体 |
| `intracellular.bacillus_like` | `bacillus_type` | *Bacillus* sp. CDB3 相关新型 PhaZ |

## 仓库结构

```text
PHB_gtdb/
  README.md
  RUN_MANIFEST.md
  environment.yml
  main_pipeline.sh
  scripts/
    01_phb_search.py
    01b_archaea_search.py
    02_extract_sequences.py
    03_msa.py
    04_phylogeny.py
    05_annotation.py
    06_visualization.py
    07_nature_figures.py
    08_prepare_ribosomal_hmms.py
    08_grodon_growth.py
  docs/
    PRESENTATION_SUMMARY.md
    REPRODUCIBLE_WORKFLOW.md
    SCRIPT_INDEX.md
    FIGURE_CAPTIONS.md
    GROWTH_RATE_ANALYSIS.md
    CODEX_WORKLOG.md
  figures/nature/
    figure1_workflow_funnel.pdf/png/svg
    figure2_phylum_heatmap.pdf/png/svg
    figure3_subtype_lipase.pdf/png/svg
    figure4_genera_phylogeny.pdf/png/svg
    figure5_grodon_growth_comparison_hmm_allmatched.pdf/png/svg
  figure_data/
    作图所需的轻量源数据
```

## 复现方式

在 T141 或配置了相同数据路径的 Linux 服务器上：

```bash
conda env create -f environment.yml
conda activate phb_gtdb

export PHB_PROJECT_ROOT=/home/data/haoyu/PHB_gtdb
export GTDB_ROOT=/home/data/haoyu/GTDB

bash main_pipeline.sh --threads 30
python scripts/check_results.py
python scripts/07_nature_figures.py
```

gRodon2 生长速率分析：

```bash
Rscript scripts/install_grodon2.R /home/data/haoyu/software/gRodon2
python scripts/08_prepare_ribosomal_hmms.py

nohup python scripts/08_grodon_growth.py \
  --heg-method hmm \
  --threads 8 \
  --marker-threads 1 \
  --max-per-genus 0 \
  --resume \
  --output-label hmm_allmatched \
  > results/logs/grodon_growth_hmm_allmatched.out 2>&1 &
```

实时可视化监控：

```bash
python scripts/09_monitor_grodon_progress.py \
  --interval 60 \
  --history results/tables/grodon_growth_monitor_history.tsv \
  --plot results/figures/grodon_growth_progress.png
```

详细步骤、每个脚本的输入输出和原理见：

- [docs/REPRODUCIBLE_WORKFLOW.md](docs/REPRODUCIBLE_WORKFLOW.md)
- [docs/SCRIPT_INDEX.md](docs/SCRIPT_INDEX.md)
- [docs/GROWTH_RATE_ANALYSIS.md](docs/GROWTH_RATE_ANALYSIS.md)
- [RUN_MANIFEST.md](RUN_MANIFEST.md)

## 已生成论文图

| Figure | 文件 | 主题 |
|---|---|---|
| Figure 1 | [figure1_workflow_funnel.pdf](figures/nature/figure1_workflow_funnel.pdf) | 整体分析流程与数据筛选漏斗 |
| Figure 2 | [figure2_phylum_heatmap.pdf](figures/nature/figure2_phylum_heatmap.pdf) | 门水平 PhaZ 分布与亚型热图 |
| Figure 3 | [figure3_subtype_lipase.pdf](figures/nature/figure3_subtype_lipase.pdf) | 五类 PhaZ 亚型组成与 lipase box 验证 |
| Figure 4 | [figure4_genera_phylogeny.pdf](figures/nature/figure4_genera_phylogeny.pdf) | Top 属分布与系统发育树概览 |
| Figure 5 | [figure5_grodon_growth_comparison_hmm_allmatched.pdf](figures/nature/figure5_grodon_growth_comparison_hmm_allmatched.pdf) | gRodon2 同属 phaZ+ vs phaZ- 预测最大生长速率比较 |

图题与图注见 [docs/FIGURE_CAPTIONS.md](docs/FIGURE_CAPTIONS.md)。

## 结果文档

- [docs/PRESENTATION_SUMMARY.md](docs/PRESENTATION_SUMMARY.md): 汇报用完整总结
- [docs/RESULTS.md](docs/RESULTS.md): 统计结果汇总
- [docs/METHODS.md](docs/METHODS.md): 方法细节
- [docs/PIPELINE.md](docs/PIPELINE.md): 原始分析流程说明
- [docs/PHAZ_REFERENCES.md](docs/PHAZ_REFERENCES.md): 14 条 PhaZ 参考序列
- [docs/REVIEW.md](docs/REVIEW.md): 项目检查与改进建议

## 数据追踪原则

- `data/raw/`, `data/processed/`, `data/external/` 不上传 GitHub，因为包含 GTDB 原始数据、外部数据库和可再生中间文件。
- `results/tables/` 和 `results/logs/` 默认不上传，因为服务器可重新生成，且可能持续更新。
- `figures/nature/` 和 `figure_data/` 已纳入仓库，用于汇报和论文图复现。
- 所有关键运行参数、服务器路径和预期结果数量记录在 [RUN_MANIFEST.md](RUN_MANIFEST.md)。
