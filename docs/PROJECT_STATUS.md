---
name: project-complete-status
description: "Complete project snapshot — all results, scripts, figures, pipeline steps, and task status as of 2026-07-02"
metadata: 
  node_type: memory
  type: project
  originSessionId: 040b394a-9136-44a8-8771-be79b2ce420d
---

# PHB_gtdb 项目完整状态快照

**更新日期**: 2026-07-02 | **总提交数**: ~25 | **GitHub**: https://github.com/anzhiw2-gif/PHB_gtdb

---

## 1. 核心结果摘要

| 指标 | 数值 |
|---|---|
| 搜索基因组 | 165,468 细菌 + 10,122 古菌 |
| DIAMOND 原始命中 | 7,068 基因组, 16,486 hits |
| CD-HIT c95 去冗余 | 8,731 序列, 7,051 基因组 |
| 长度过滤 ≥100aa | 7,478 序列, 6,224 基因组 |
| **验证集（论文核心）** | **6,532 序列, 5,531 基因组** |
| 覆盖 | 13 门, 1,135 属 |
| 古菌 PhaZ | 0（确认无） |
| 平均拷贝数 | 1.18/基因组（84.7% 单拷贝） |

## 2. PhaZ 五亚型

| 亚型 | 验证序列 | 验证基因组 | Lipase box | 信号肽 |
|---|---|---|---|---|
| intracellular.cupriavidus_like | 3,668 | 3,358 | 10.1% | 7.8% |
| intracellular.ralstonia_like | 2,221 | 2,062 | 9.8% | 8.6% |
| intracellular.bacillus_like | 5 | 5 | 50.0% | 0% |
| extracellular.general | 404 | 365 | 28.6% | 7.9% |
| extracellular.lemoignei_like | 234 | 232 | 53.2% | 16.7% |

## 3. gRodon2 生长速率

| 指标 | 值 |
|---|---|
| 分析基因组 | 8,788 (4,394+ / 4,394-) |
| 覆盖属 | 899 |
| 属均值 δ | -0.0027/h (95% CI: -0.0155, +0.0093) |
| Wilcoxon p | **0.459（不显著）** |
| 效应量 r | 0.025（可忽略） |
| 结论 | **同属内 phaZ+ vs phaZ- 无显著生长速率差异** |

## 4. 补充分析（Chapter 2 论文用）

| # | 任务 | 状态 | 输出 |
|---|---|---|---|
| T1 | 拷贝数分布 | ✅ | 84.7% 单拷贝，验证集 vs c95 对比 |
| T2 | PhaC-phaZ 共现 | ✅ | 12 属同时含 phaZ+phaC |
| T3 | 3D 结构预测 | ✅ | ESM Atlas API，每亚型 5 条 PDB |
| T4 | 催化域保守性 | ✅ | 5 亚型位置保守谱 CSV |
| T5 | 信号肽预测 | ✅ | Kyte-Doolittle 疏水性 |
| T6 | 宏基因组种子 | ✅ | 45 条代表序列 FASTA |
| T7 | 汇总出图 | 🔄 | Fig 2.6-2.8 待生成 |

## 5. 脚本索引

### 核心管道 (01-06)
| 脚本 | 功能 |
|---|---|
| `01_phb_search.py` | 细菌 PhaZ 搜索（Pyrodigal + DIAMOND） |
| `01b_archaea_search.py` | 古菌 PhaZ 搜索 |
| `02_extract_sequences.py` | 提取蛋白序列 + CD-HIT + 过滤 |
| `03_msa.py` | 亚型分型 + MAFFT + trimAl |
| `04_phylogeny.py` | IQ-TREE/FastTree 建树 |
| `05_annotation.py` | 功能注释（Pfam, eggNOG, KOfam） |
| `06_visualization.py` | 常规可视化 |

### 论文图 (07)
| 脚本 | 功能 |
|---|---|
| `07_nature_figures.py` | Nature 风格 Fig 1-5 生成 |

### gRodon2 扩展 (08-09)
| 脚本 | 功能 |
|---|---|
| `08_grodon_growth.py` | gRodon2 批量分析主脚本 |
| `08_make_grodon_seed.R` | 核糖体 seed 生成 |
| `08_prepare_ribosomal_hmms.py` | Pfam 核糖体 HMM 准备 |
| `08_run_grodon_one.R` | 单基因组 gRodon2 包装 |
| `09_monitor_grodon_progress.py` | 实时监控 |
| `10_balance_grodon_by_genus.py` | 属内配对平衡 |
| `11_grodon_growth_stats.py` | 统计检验 |

### Chapter 2 补充分析 (t1-t6)
| 脚本 | 任务 | 功能 |
|---|---|---|
| `t1_validated.py` | T1 | 验证集拷贝数分布 |
| `t2_phac_search.py` | T2 | PhaC-phaZ 共现分析 |
| `t3_select_reps.py` | T3 | 代表序列选拔 |
| `t3_esm_api.py` | T3 | ESM Atlas API 3D 折叠 |
| `t3_esm2_contacts.py` | T3 | ESM-2 接触图预测 |
| `t4_conservation.py` | T4 | 催化域保守性 |
| `t5_signal_peptide.py` | T5 | Kyte-Doolittle 信号肽 |
| `t6_seed_sequences.py` | T6 | 宏基因组种子序列 |

### 工具
| 脚本 | 功能 |
|---|---|
| `config.py` | 全局路径配置 |
| `utils.py` | FASTA/日志/命令工具函数 |
| `check_results.py` | 结果完整性检查 |

## 6. 图件清单

| 图 | 文件 | 内容 | 状态 |
|---|---|---|---|
| Fig 1 | figure1_workflow_funnel | 分析流程 + 筛选漏斗 | ✅ |
| Fig 2 | figure2_phylum_heatmap | 门分布 + 亚型热图 | ✅ |
| Fig 3 | figure3_subtype_lipase | 亚型组成 + Lipase box | ✅ |
| Fig 4 | figure4_genera_phylogeny | 属分布 + 系统发育树 | ✅ |
| Fig 5 | figure5_grodon_growth_comparison | gRodon2 生长速率 | ✅ |
| Fig 2.6 | — | 拷贝数 + PhaC 共现 | 🔄 待生成 |
| Fig 2.7 | — | 3D 结构 + 催化域 | 🔄 待生成 |
| Fig 2.8 | — | 信号肽 + 种子序列 | 🔄 待生成 |

## 7. 文档清单

| 文档 | 内容 |
|---|---|
| `README.md` | 项目主文档 |
| `RUN_MANIFEST.md` | 运行清单与预期结果 |
| `GENOME_VS_PROTEIN_COUNTS.md` | 基因组/序列数量完整对比 |
| `PRESENTATION_SUMMARY.md` | 汇报用总结 |
| `FINAL_REPORT.md` | 最终报告 |
| `METHODS.md` | 方法细节（论文用） |
| `RESULTS.md` | 统计结果汇总 |
| `PIPELINE.md` | 分析流程说明 |
| `PHAZ_REFERENCES.md` | 14 条参考序列 |
| `SCRIPT_INDEX.md` | 脚本索引 |
| `FIGURE_CAPTIONS.md` | 图题图注（5 张） |
| `REPRODUCIBLE_WORKFLOW.md` | 可复现流程 |
| `GROWTH_RATE_ANALYSIS.md` | gRodon2 设计文档 |
| `CODEX_WORKLOG.md` | Codex 协作记录 |
| `REVIEW.md` | 管线评审 |
| `SEARCH_REPORT.md` | 搜索结果报告 |
| `superpowers/specs/2026-07-02-chapter2-supplement-design.md` | Chapter 2 设计文档 |

## 8. 运行环境

| 项目 | 值 |
|---|---|
| 主要服务器 | **T141** (10.16.1.141) |
| 项目路径 | `/home/data/haoyu/PHB_gtdb` |
| GTDB 数据 | `/home/data/haoyu/GTDB/` |
| Python 环境 | `/home/data/haoyu/miniconda3/envs/phb_gtdb` |
| R/gRodon2 环境 | `/home/data/haoyu/miniconda3/envs/grodon2` |
| ESMFold 环境 | `/home/data/haoyu/miniconda3/envs/esmfold` |
| ColabFold 环境 | `/home/data/haoyu/miniconda3/envs/colabfold` |
| 本地开发 | `D:/PHB_gtdb` (Windows) |
| GTDB 版本 | R232 |
| CPU | 80 核 |
| 内存 | 1 TB |
| GPU | ❌ 无 |

## 9. 关键数据文件路径 (T141)

```
data/processed/
├── phb_search_results.tsv           # 7,068 行 DIAMOND 搜索结果
├── phaz_proteins_all.fasta          # 8,767 条提取序列
├── phaz_proteins_c95.fasta          # 8,731 条 c95 去冗余
├── phaz_proteins_filtered.fasta     # 7,478 条 ≥100aa
├── phaz_proteins_validated.fasta    # 6,532 条验证集 ★
├── phaz_intracellular.fasta         # 4,693 条 (c95)
├── phaz_ralstonia.fasta             # 2,948 条 (c95)
├── phaz_bacillus_type.fasta         # 60 条 (c95)
├── phaz_extracellular.fasta         # 596 条 (c95)
├── phaz_extracellular_lemoignei.fasta # 434 条 (c95)
├── phaz_*_trim.fasta                # 8,178 条 trimAl 修剪后
├── phaz_3d_representatives.fasta    # 5 条 3D 代表序列
└── phaz_metagenome_seeds.fasta      # 45 条宏基因组种子

results/
├── tables/
│   ├── grodon_growth_*.tsv          # gRodon2 全量结果
│   ├── phaz_validated_*.tsv         # 验证集统计表
│   ├── phaz_signal_peptide_*.tsv    # 信号肽预测
│   ├── phaz_phac_cross_reference.tsv # PhaC 共现
│   ├── phaz_metagenome_seeds.tsv    # 种子序列表
│   └── conservation_*.tsv           # 催化域保守性
├── structures/
│   ├── intra_cupriavidus_*.pdb      # 胞内 Cupriavidus 型 3D
│   ├── intra_ralstonia_*.pdb        # 胞内 Ralstonia 型 3D
│   ├── intra_bacillus_*.pdb         # 胞内 Bacillus 型 3D
│   ├── extra_general_*.pdb          # 胞外一般型 3D
│   └── extra_lemoignei_*.pdb        # 胞外 Lemoignei 型 3D
├── tree/                            # 5 亚型系统发育树
└── figures/                         # 常规可视化图
```

## 10. 论文使用指南

| 场合 | 数字 |
|---|---|
| PhaZ 蛋白总数 | **6,532** |
| 含 phaZ 的基因组 | **5,531** |
| 分布门/属 | **13 门, 1,135 属** |
| 平均拷贝数 | **1.18** |
| 系统发育建树 | c95 级 8,178 条 |
| gRodon2 分析 | 8,788 基因组, 899 属 |
| 3D 结构 | 25 条 (每亚型 5 条) |

[[project-overview]] [[key-results]] [[technical-stack]] [[server-config]] [[phaz-subtypes]] [[methodology-decisions]] [[grodon2-results]] [[user-profile]]
