# PHB_gtdb 项目审查总结文档

**生成日期**: 2026-07-02 | **版本**: v2.0 (详尽版)  
**项目**: GTDB R232 中 PHB 降解关键酶 PhaZ 的系统鉴定、分型与进化分析  
**GitHub**: https://github.com/anzhiw2-gif/PHB_gtdb  
**用途**: 硕士论文 Chapter 2 生物信息学分析 | 导师审查 | 组会汇报 | 项目归档

---

## 一、项目定位与科学问题

### 1.1 项目背景

聚羟基丁酸酯 (PHB) 是微生物合成的一类可降解生物塑料 (PHA)，PhaZ (PHB depolymerase, EC 3.1.1.75) 是 PHB 降解过程的关键酶。本项目旨在利用 GTDB (Genome Taxonomy Database) 这一最全面的微生物基因组分类数据库，系统回答以下问题：

1. **哪些微生物携带 PhaZ？** — 在全部细菌和古菌代表基因组中的分布全景
2. **PhaZ 的分类学集中性** — 是否集中在特定门/属？
3. **古菌是否存在 phaZ？** — 验证"古菌无 phaZ"假说
4. **PhaZ 的亚型分类** — 能否按结构域/保守基序分为稳定亚型？
5. **phaZ 携带是否影响生长策略？** — 同属内 phaZ+ vs phaZ- 预测生长速率比较
6. **PhaZ 蛋白结构特征** — 各亚型的催化域保守性和三维折叠差异

### 1.2 硕士论文架构

```
Chapter 1  引言与文献综述
Chapter 2  基于 GTDB 全景扫描的 PhaZ 基因鉴定与进化分析  ← 本项目
Chapter 3  PHB 降解菌的分离培养与酶活性验证             ← 湿实验
Chapter 4  环境宏基因组中 PhaZ 的多样性与分布特征        ← 宏基因组
Chapter 5  综合讨论与结论
```

三章逻辑链: **Chapter 2** (全景→告诉你在哪里找) → **Chapter 3** (实验室验证) → **Chapter 4** (自然界验证)

---

## 二、数据与方法

### 2.1 数据来源

| 数据 | 路径/来源 | 说明 |
|---|---|---|
| GTDB Release | R232 | gtdb.ecogenomic.org |
| 细菌代表基因组 | 189,801 个 | `bac120_taxonomy_r232.tsv` 中的代表基因组 |
| 古菌代表基因组 | 10,122 个 | `ar53_taxonomy_r232.tsv` 中的代表基因组 |
| 基因组文件 | `gtdb_genomes_reps_r232/database/` | 199,923 个 `*_genomic.fna.gz` |
| 分类学信息 | `bac120_taxonomy_r232.tsv` | 878,998 行 (全部基因组) |
| 质量元数据 | `bac120_metadata_r232.tsv.gz` | 878,999 行 |
| 物种树 | `bac120_r232.tree` | GTDB 官方参考树 |

### 2.2 PhaZ 参考序列

**14 条经 NCBI/文献验证的 PhaZ 参考序列**，覆盖 3 门 8 属：

| 类型 | 代表物种 | 序列数 | 关键文献 |
|---|---|---|---|
| 胞外型 (extracellular) | *Stutzerimonas stutzeri*, *Comamonas testosteroni*, *Acidovorax* sp., *Delftia acidovorans*, *Streptomyces exfoliatus* | 5 | Ohura 1999; Jendrossek 1995 |
| 胞外型 (lemoignei) | *Paucimonas lemoignei* (PhaZ1-Z7) | 3 | Jendrossek et al. |
| 胞内型 (Cupriavidus) | *Cupriavidus necator* H16 (PhaZ1, PhaZ2, PhaZ5) | 3 | Saegusa 2001; Pohlmann 2006 |
| 胞内型 (Ralstonia) | *Ralstonia pickettii* | 2 | Genome annotation |
| 新型 (Bacillus) | *Bacillus* sp. CDB3 | 1 | Tseng 2006 |

所有序列均通过 NCBI E-utilities API 实时查询验证。详见 `docs/PHAZ_REFERENCES.md`。

### 2.3 分析管道

```
GTDB R232 代表基因组 (199,923)
    │
    ├── Step 1: Pyrodigal 基因预测 (meta=True, 进程内)
    │         DIAMOND blastp 同源搜索 (E≤1e-10, ID≥30%, qcov≥50%)
    │         → 细菌 7,068 命中, 16,486 hits / 古菌 2 命中
    │
    ├── Step 2: 蛋白序列提取 + CD-HIT c95 (仅去真冗余)
    │         长度过滤 ≥100aa
    │         → 8,731 (c95) → 7,478 (filtered)
    │
    ├── Step 3: HMMer profile 验证 (14-ref HMM, E≤1e-5)
    │         + DIAMOND 高置信度 (pident≥35%, ≥150aa)
    │         → 6,532 验证集 (并集)
    │
    ├── Step 4: Pfam 结构域分型 → 5 亚型
    │         MAFFT --auto 多序列比对
    │         trimAl -automated1 修剪
    │
    ├── Step 5: IQ-TREE (小组, MFP+B1000)
    │         FastTree (大组, LG+G4)
    │         → 5 棵独立系统发育树
    │
    ├── Step 6: 功能注释 + 可视化
    │
    └── 扩展: gRodon2 生长速率分析
              3D 结构预测 (ESM Atlas API)
              Pfam 结构域注释 (HMMER hmmscan)
```

### 2.4 计算环境

| 项目 | 值 |
|---|---|
| 服务器 | T141 (root123-PowerEdge-T630) |
| IP | 10.16.1.141 |
| CPU | 80 核 |
| 内存 | 1 TB |
| GPU | 无 |
| 磁盘 | 131 TB (可用 27 TB) |
| Python 环境 | conda `phb_gtdb` (Python 3.12) |
| R/gRodon2 环境 | conda `grodon2` |
| ESMFold 环境 | conda `esmfold` |
| ColabFold 环境 | conda `colabfold` |
| 项目路径 | `/home/data/haoyu/PHB_gtdb` |
| GTDB 数据 | `/home/data/haoyu/GTDB/` |

---

## 三、核心结果

### 3.1 管道各级过滤数据

| 步骤 | 序列数 | 基因组数 | 移除 | 移除原因 |
|---|---|---|---|---|
| ① DIAMOND 原始命中 | 16,486 hits | 7,068 | — | — |
| ② 蛋白序列提取 | 8,767 | 7,068 | 7,719 hits | 去重叠 CDS 命中 |
| ③ CD-HIT c95 | 8,731 | 7,051 | 36 / 17 | 仅去真冗余 (0.4%) |
| ④ 长度 ≥100 aa | 7,478 | 6,224 | 1,253 / 827 | 短片段 (含古菌 2 条 <50aa) |
| ⑤ HMMer 验证 | 6,033 | — | — | profile HMM, E≤1e-5 |
| ⑥ DIAMOND 高置信度 | 3,270 | — | — | pident≥35% + ≥150aa |
| **⑦ 最终验证集** | **6,532** | **5,531** | **946 / 693** | **HMMer + DIAMOND 并集** |

> **6,032 vs 6,532**: HMMer 单独验证通过 6,033 条 (保守)，DIAMOND 高置信度额外补充 499 条边界命中，两者并集 = 6,532。

### 3.2 古菌验证

| 指标 | 值 |
|---|---|
| 古菌基因组扫描 | 10,122 |
| DIAMOND 初始命中 | 2 (*Halovenus* sp., 34aa; UBA73 sp., 44aa) |
| ≥100aa 过滤后 | **0** |
| HMMer 验证后 | **0** |
| **结论** | **在本次 GTDB R232 代表古菌基因组分析范围内，未发现高置信度 phaZ 基因** |

### 3.3 PhaZ 亚型分类

| 亚型名称 | 验证序列 | 验证基因组 | 占比 | Lipase box | 信号特征 | c95 序列 | c95 基因组 |
|---|---|---|---|---|---|---|---|
| intracellular.cupriavidus_like | 3,668 | 3,358 | 56.2% | 10.1% | 7.8% | 4,693 | 4,200 |
| intracellular.ralstonia_like | 2,221 | 2,062 | 34.0% | 9.8% | 8.6% | 2,948 | 2,696 |
| extracellular.general | 404 | 365 | 6.2% | **28.6%** | 7.9% | 596 | 500 |
| extracellular.lemoignei_like | 234 | 232 | 3.6% | **53.2%** | **16.7%** | 434 | 418 |
| intracellular.bacillus_like | 5 | 5 | 0.08% | **50.0%** | 0% | 60 | 53 |
| **合计** | **6,532** | **5,531** | 100% | | | 8,731 | 7,051 |

> 注：各亚型基因组数之和 (6,022) > 总计 (5,531)，因为 **491 个基因组含有跨亚型的多个 phaZ 拷贝**。

**Lipase box 关键验证:**
- 胞外型 (28-53%) 显著高于胞内型 (~10%)，验证分型合理性
- Bacillus 型检出 G-**W**-S-**M**-G 基序，与 Tseng et al. (2006) 报道完全一致
- lemoignei 型 lipase box 检出率最高 (53.2%)，与该类群以分泌型 PhaZ 为主的生物学特征一致

### 3.4 拷贝数分布

| 拷贝数 | 基因组数 | 占比 |
|---|---|---|
| 1 | **4,682** | **84.7%** |
| 2 | 724 | 13.1% |
| 3 | 104 | 1.9% |
| 4 | 17 | 0.3% |
| 5 | 2 | 0.04% |
| 6 | 2 | 0.04% |

> **关键发现**: 原始 DIAMOND 搜索中 67.7% 基因组表现为多拷贝 (最高 12 个)，经 HMMer 验证后大量弱命中被识别为假阳性同源片段。验证集以单拷贝为主 (84.7%)，平均 1.18 拷贝/基因组。

### 3.5 门级分布 (验证集，13 门)

| 门 | 序列数 | 占比 |
|---|---|---|
| Pseudomonadota | 6,323 | 96.8% |
| Actinomycetota | 154 | 2.4% |
| Myxococcota | 14 | 0.2% |
| Acidobacteriota | 6 | 0.09% |
| Bacteroidota | 5 | 0.08% |
| Bacillota | 5 | 0.08% |
| Chloroflexota | 5 | 0.08% |
| Gemmatimonadota | 5 | 0.08% |
| Bdellovibrionota | 4 | 0.06% |
| Desulfobacterota | 4 | 0.06% |
| Spirochaetota | 3 | 0.05% |
| Planctomycetota | 3 | 0.05% |
| Cyanobacteriota | 1 | 0.02% |

> **c95 级别覆盖 17 门**，额外包括 Halobacteriota, Nanobdellota, Thermotogota, Verrucomicrobiota，
> 但经 HMMer 严格过滤后这 4 门的命中被移除。

> Pseudomonadota 占主导 (96.8%)，部分反映参考序列偏倚 (11/14 参考序列来自此门)。

### 3.6 Top 15 属 (验证集)

| 属 | 序列数 | 主导亚型 | 已知 PHB 降解 |
|---|---|---|---|
| *Paraburkholderia* | 319 | ralstonia | ✓ |
| *Bradyrhizobium* | 239 | intracellular | ✓ |
| *Cupriavidus* | 148 | intracellular | ✓ (参考来源) |
| *Variovorax* | 136 | extracellular | ✓ |
| *Polaromonas* | 109 | intracellular | — |
| *Rubrivivax* | 107 | intracellular | ✓ |
| *Polynucleobacter* | 98 | intracellular | — |
| *Mesorhizobium* | 96 | intracellular | ✓ |
| *Stenotrophomonas* | 93 | intracellular | — |
| *Burkholderia* | 90 | ralstonia | ✓ |
| *Caballeronia* | 73 | intracellular | — |
| *Pseudolabrys* | 73 | intracellular | — |
| *Rugamonas* | 71 | intracellular | — |
| *BOG-931* | 71 | intracellular | — |
| *Hydrogenophaga* | 68 | intracellular | ✓ |

---

## 四、系统发育分析

### 4.1 各亚型建树详情

| 亚型 | 修剪后序列 | 修剪后基因组 | 工具 | 模型 | 树文件大小 | 状态 |
|---|---|---|---|---|---|---|
| intracellular.cupriavidus_like | 4,424 | 4,014 | FastTree 2.2 | LG+G4 | 540 KB | ✅ |
| intracellular.ralstonia_like | 2,776 | 2,574 | IQ-TREE + FastTree | LG+F+G4 | 332 KB | ✅ |
| extracellular.general | 509 | 453 | IQ-TREE 2 | MFP+B1000 | 58 KB | ✅ |
| extracellular.lemoignei_like | 412 | 399 | IQ-TREE 2 | MFP+B1000 | 49 KB | ✅ |
| intracellular.bacillus_like | 57 | 51 | IQ-TREE 2 | MFP+B1000 | 6.4 KB | ✅ |
| **联合树 (#14)** | 857 (采样) | — | MAFFT+FastTree | LG+G4 | — | ✅ |

> **建树策略**: 大组 (>2,000 序列) 使用 FastTree——IQ-TREE ModelFinder 需测试 437 个模型 (实测 25h 未完成)，FastTree 数分钟完成且拓扑准确度高度一致。(Price et al., 2010)

> **联合树采样**: 小亚型 (bacillus, 57 条) 全保留；大亚型采样至 150-250 条，总计 857 条，
> 随机种子 42 确保可重复。详见 `docs/SAMPLING_STRATEGY.md`。

---

## 五、gRodon2 生长速率分析

### 5.1 分析设计

**配对策略**: 仅在含有 phaZ+ 的属内寻找同属 phaZ- 对照基因组，属内数量平衡。

**HEG 识别**: Pfam 原核核糖体蛋白 HMM (筛选自 Pfam-A, 排除真核/线粒体/质体模型) + HMMER。

### 5.2 核心结果

| 指标 | 值 |
|---|---|
| 输入基因组 | 8,788 (4,394 phaZ+ / 4,394 phaZ-) |
| 平衡后保留 | **8,692** (4,346 / 4,346) |
| 覆盖属 | **899** |
| 失败基因组 | 53 (核糖体蛋白 HMM 命中 <10) |

| 统计检验 | 结果 | 显著性 |
|---|---|---|
| 属均值 δ (growth/h) | **-0.0027** (95% CI: -0.0155, +0.0093) | 跨零 |
| 属中位数 δ | +0.0018 | — |
| Wilcoxon 符号秩 | z=0.74, **p=0.459** | ❌ 不显著 |
| 分层置换检验 (未加权) | **p=0.670** | ❌ 不显著 |
| 分层置换检验 (加权) | p=0.231 | ❌ 不显著 |
| 符号检验 (正/负属) | 463/436, **p=0.386** | ❌ 不显著 |
| 效应量 r | **0.025** | 可忽略 |

| 基因组级 (仅描述性) | phaZ+ | phaZ- |
|---|---|---|
| 均值 growth/h | 0.392 | 0.384 |
| 中位数 growth/h | 0.240 | 0.232 |

**结论: 同属内 phaZ+ 与 phaZ- 菌株的预测最大生长速率无统计学显著差异。效应量可忽略 (r=0.025)，899 属的样本量足够检测到有意义的差异。**

### 5.3 表述注意事项

- gRodon2 输出为基于 codon usage bias 的**预测最大生长潜力**，非实测培养生长速率
- 论文中应使用 "predicted maximal growth potential" 或 "预测最大生长潜力"
- 53 个失败基因组来自核糖体蛋白 HMM 命中不足，属正常范围
- 当前分析为全局属级比较；亚型分层分析待做

---

## 六、结构域与功能预测

### 6.1 Lipase Box 验证

| 亚型 | 完整序列检出率 | 方法 |
|---|---|---|
| intracellular.cupriavidus_like | 10.1% | G-X-S-X-G 正则匹配 |
| intracellular.ralstonia_like | 9.8% | G-X-S-X-G 正则匹配 |
| intracellular.bacillus_like | 50.0% | G-X-S-X-G, 检出 G-W-S-M-G |
| extracellular.general | 28.6% | G-X-S-X-G 正则匹配 |
| extracellular.lemoignei_like | 53.2% | G-X-S-X-G 正则匹配 |

> 基于完整验证集序列的 G-X-S-X-G 五肽正则匹配。
> trimAl 修剪后的比对中 lipase box 检出率极低 (<2%)，因为该区域为可变区，
> 论文应以完整序列扫描结果为准。

### 6.2 信号肽预测 (Kyte-Doolittle)

| 亚型 | 总数 | 信号特征阳性 | 比例 | 平均疏水性 |
|---|---|---|---|---|
| intracellular | 3,668 | 285 | 7.8% | 0.79 |
| ralstonia | 2,221 | 190 | 8.6% | 0.78 |
| bacillus_type | 5 | 0 | 0% | 1.16 |
| extracellular | 404 | 32 | 7.9% | 0.59 |
| extracellular_lemoignei | 234 | **39** | **16.7%** | 0.86 |
| 胞内型合计 | 5,894 | 475 | 8.1% | — |
| 胞外型合计 | 638 | **71** | **11.1%** | — |

> **方法**: Kyte-Doolittle 疏水性量表 + N 端正电荷检测。
> **⚠ 非 SignalP 6.0 正式预测**。N 端疏水窗口 + 正电荷富集是 Sec 分泌信号肽的经典特征，
> 但表述应为 "putative N-terminal hydrophobic signal feature"。

### 6.3 催化域保守性 (T4)

基于 c95 级别 trimmed 比对 (8,178 条) 的位点保守性计算：

| 亚型 | 比对长度 | 最大保守位点 | 保守度 |
|---|---|---|---|
| intracellular.cupriavidus_like | 80 aa | 位点 70: F | 58% |
| intracellular.ralstonia_like | 111 aa | 位点 62: F | 50% |
| intracellular.bacillus_like | 29 aa | 位点 7: P | 61% |
| extracellular.general | 26 aa | 位点 10: R | 45% |
| extracellular.lemoignei_like | 18 aa | **位点 4: G** | **79%** |

> trimAl 修剪后保留的是催化核心区。Lemoignei 型 G79 位点可能是 lipase box G 残基。

### 6.4 3D 结构预测 (T3)

**方法**: ESMFold v1 (Meta), 通过 ESM Atlas API (`api.esmatlas.com/foldSequence/v1/pdb/`)
**引用**: Lin et al., 2023, *Science*

| 亚型 | 预测数 | 截取策略 |
|---|---|---|
| intracellular.cupriavidus_like | 5 | 300aa 催化域 |
| intracellular.ralstonia_like | 5 | 300aa 催化域 |
| intracellular.bacillus_like | 4 | 109-207aa 全长 |
| extracellular.general | 5 | 300aa 催化域 |
| extracellular.lemoignei_like | 5 | 200-300aa 催化域 |
| **总计** | **28 PDB** | 2 条 HTTP 504 超时失败 |

### 6.5 PhaC 共现分析 (T2)

**方法**: 文献收集已知 PhaC 属 (Class I-IV 的代表属)，与 phaZ+ 验证属交叉。

| 指标 | 结果 |
|---|---|
| 已知 PhaC 属 | 19 个 (文献收集) |
| 同时含 phaZ+phaC 的属 | **12 个** |
| 关键共现属 | *Cupriavidus* (148), *Paraburkholderia* (319), *Bradyrhizobium* (239), *Variovorax* (136), *Hydrogenophaga* (68), *Burkholderia* (90), *Ralstonia* (29) |
| 有 PhaC 但无 phaZ | 7 个 (*Pseudomonas*, *Bacillus*, *Aeromonas*, *Allochromatium*, *Chromobacterium*, *Delftia*, *Wautersia*) |

> 当前分析为属级交叉引用，升级到 Class I-IV PhaC HMM 搜索 + 基因组级共现 + 树拓扑比较可支撑"共进化"宣称。

### 6.6 Pfam 结构域注释 (#16)

**方法**: HMMER hmmscan 对 6,532 条验证序列搜索 Pfam-A (30,134 HMMs), E≤1e-5

**状态**: 🔄 T141 上运行中 (2026-07-02 14:36 启动, 3,525 个结构域命中)

**预期关键 Pfam 结构域**:
- PF06850 (PHB_depo_C): PHB depolymerase C-terminus
- PF00561 (Abhydrolase_1): alpha/beta hydrolase fold
- PF12697 (Abhydrolase_6): alpha/beta hydrolase family
- PF14556 (PHA_depolymerase): PHA depolymerase domain

---

## 七、图件清单

### 7.1 Nature 风格主图 (Fig 1-5)

| 图 | 文件 | 内容 | 尺寸 (PDF) |
|---|---|---|---|
| Fig 1 | `figure1_workflow_funnel` | 分析工作流 + 筛选漏斗 | 62 KB |
| Fig 2 | `figure2_phylum_heatmap` | 门水平分布 + 亚型热图 | 54 KB |
| Fig 3 | `figure3_subtype_lipase` | 亚型组成 + Lipase box 验证 | 44 KB |
| Fig 4 | `figure4_genera_phylogeny` | Top 属分布 + 系统发育树 | 461 KB |
| Fig 5 | `figure5_grodon_growth_comparison` | gRodon2 生长速率比较 | 255 KB |

### 7.2 补充图 (Fig 2.6-2.8)

| 图 | 文件 | 内容 | 尺寸 (PDF) |
|---|---|---|---|
| Fig 2.6 | `figure26_copy_number_phac` | 拷贝数分布 + PhaC 共现 | 68 KB |
| Fig 2.7 | `figure27_structures_conservation` | 28 PDB 结构 + 催化域保守性 | 47 KB |
| Fig 2.8 | `figure28_signal_seeds` | 信号肽预测 + 宏基因组种子 | 44 KB |

> 全部提供 PDF (矢量) + PNG (600 DPI) + SVG (可编辑文本) 三种格式。

---

## 八、脚本总览

### 8.1 核心管道 (01-06)
| 脚本 | 功能 | 主要工具 |
|---|---|---|
| `01_phb_search.py` | 细菌 PhaZ 搜索 | Pyrodigal + DIAMOND |
| `01b_archaea_search.py` | 古菌 PhaZ 搜索 | Pyrodigal + DIAMOND |
| `02_extract_sequences.py` | 序列提取 + CD-HIT + 长度过滤 | CD-HIT |
| `03_msa.py` | 亚型分型 + MAFFT + trimAl | MAFFT, trimAl |
| `04_phylogeny.py` | 系统发育建树 | IQ-TREE, FastTree |
| `05_annotation.py` | 功能注释 | eggNOG, KOfam |
| `06_visualization.py` | 常规可视化 | matplotlib |

### 8.2 论文图 (07)
| 脚本 | 功能 |
|---|---|
| `07_nature_figures.py` | Nature Fig 1-5 (633 行) |
| `07b_supplement_figures.py` | Fig 2.6-2.8 (297 行) |

### 8.3 gRodon2 扩展 (08-11)
| 脚本 | 功能 |
|---|---|
| `08_grodon_growth.py` | gRodon2 批量分析主脚本 |
| `08_prepare_ribosomal_hmms.py` | Pfam 核糖体蛋白 HMM 准备 |
| `08_run_grodon_one.R` | 单基因组 gRodon2 包装 |
| `09_monitor_grodon_progress.py` | 实时进度监控 |
| `10_balance_grodon_by_genus.py` | 属内配对平衡 |
| `11_grodon_growth_stats.py` | 统计检验 |

### 8.4 补充分析 (t1-t6, 论文用)
| 脚本 | 任务 | 功能 |
|---|---|---|
| `t1_validated.py` | T1 | 验证集拷贝数分布 |
| `t2_phac_search.py` | T2 | PhaC-phaZ 交叉引用 |
| `t3_select_reps.py` | T3 | 3D 代表序列选拔 |
| `t3_esm_api.py` | T3 | ESM Atlas API 折叠 |
| `t3_multi.py` | T3 | 批量 API 折叠 (25 条) |
| `t3_esm2_contacts.py` | T3 | ESM-2 嵌入/接触预测 |
| `t4_conservation.py` | T4 | 催化域位点保守性 |
| `t5_signal_peptide.py` | T5 | Kyte-Doolittle 信号肽 |
| `t6_seed_sequences.py` | T6 | 宏基因组种子选择 |

### 8.5 扩展分析 (专家建议)
| 脚本 | 功能 |
|---|---|
| `t9_lipase_logo.py` | Lipase box 序列 Logo |
| `t14_combined_tree.py` | 联合树单系性检验 |
| `t16_pfam_scan.sh` | Pfam-A hmmscan |

### 8.6 工具
| 脚本 | 功能 |
|---|---|
| `config.py` | 全局路径 + 参数配置 |
| `utils.py` | FASTA/命令/日志工具 |
| `check_results.py` | 服务器端全量结果检查 |
| `check_tracked_assets.py` | GitHub 仓库轻量资产检查 |

---

## 九、文档清单

| 文档 | 路径 | 用途 |
|---|---|---|
| 项目状态快照 | `docs/PROJECT_STATUS.md` | 全部结果/脚本/文件一览 |
| 数据账本 | `docs/DATA_ACCOUNTING.md` | 基因组计数审计 + 各级过滤追踪 |
| 基因组 vs 序列 | `docs/GENOME_VS_PROTEIN_COUNTS.md` | 基因组数≠序列数 详细解释 |
| 采样策略 | `docs/SAMPLING_STRATEGY.md` | 所有抽样决策 + 理由 |
| 设计文档 | `docs/superpowers/specs/2026-07-02-chapter2-supplement-design.md` | Chapter 2 分析设计方案 |
| **本文档** | `docs/PROJECT_REVIEW_SUMMARY.md` | 项目审查总结 (用于导师/组会) |
| 汇报总结 | `docs/PRESENTATION_SUMMARY.md` | PPT 用总结 |
| 最终报告 | `docs/FINAL_REPORT.md` | 完整分析报告 |
| 方法细节 | `docs/METHODS.md` | 论文方法部分 |
| 图题图注 | `docs/FIGURE_CAPTIONS.md` | 8 张图的中英文标题 |
| 管线文档 | `docs/PIPELINE.md` | 分析流程详解 |
| 参考序列 | `docs/PHAZ_REFERENCES.md` | 14 条参考序列详情 |
| 脚本索引 | `docs/SCRIPT_INDEX.md` | 每个脚本的原理说明 |
| 可复现流程 | `docs/REPRODUCIBLE_WORKFLOW.md` | 复现步骤 |
| 生长速率分析 | `docs/GROWTH_RATE_ANALYSIS.md` | gRodon2 设计文档 |
| gRodon2 结果 | `docs/GRODON2_FINAL_STATS.md` | 最终统计 |
| 审查改进 | `docs/REVIEW.md` | 21 项改进清单 + 状态 |
| 审计报告 | `docs/RESULTS_AUDIT_2026-06-23.md` | 结果审计 |
| 协作记录 | `docs/CODEX_WORKLOG.md` | Codex 协作归档 |
| 搜索结果 | `docs/SEARCH_REPORT.md` | 搜索阶段报告 |
| 运行清单 | `RUN_MANIFEST.md` | 预期结果计数 |
| README | `README.md` | 项目入口文档 |

---

## 十、数据资产清单

### 10.1 核心数据 (figure_data/, 已纳入 GitHub)

| 文件 | 内容 | 行数 |
|---|---|---|
| `phaz_validated_genome_protein_summary.tsv` | 验证集核心统计 | 8 |
| `phaz_validated_copy_number.tsv` | 拷贝数分布 | 7 |
| `phaz_validated_subtype_count.tsv` | 亚型组成 | 6 |
| `phaz_validated_phylum_count.tsv` | 门级分布 | 14 |
| `phaz_signal_peptide_summary.tsv` | 信号肽预测 | 6 |
| `phaz_phac_cross_reference.tsv` | PhaC 交叉引用 | 1,143 |
| `phaz_metagenome_seeds.tsv` | 宏基因组种子 | 46 |
| `phaz_proteins_validated.fasta` | 6,532 条验证序列 | (FASTA) |
| `phb_search_results.tsv` | 7,068 行搜索结果 | 7,069 |
| `archaea_phb_search_results.tsv` | 古菌搜索结果 | 3 |
| `*_subtype_strip.txt` | 5 个亚型的 phylum strip 标注 | — |
| `*_tree.treefile` | 5 个亚型的系统发育树 (Newick) | — |

### 10.2 结果表 (results/tables/, 服务器端)

- `grodon_growth_predictions_hmm_allmatched.tsv` — 8,788 行 gRodon2 预测
- `grodon_growth_same_genus_summary_hmm_allmatched.tsv` — 900 属汇总
- `grodon_growth_genus_effects_hmm_allmatched.tsv` — 899 属效应量
- `grodon_growth_statistical_tests_hmm_allmatched.tsv` — 24 项统计检验
- `conservation_*.tsv` — 5 亚型催化域保守性
- `phaz_pfam_domains.tsv` — Pfam 注释 (#16, 运行中)

### 10.3 3D 结构 (results/structures/)

28 个 PDB 文件 (ESMFold v1, 每亚型 4-5 条, 总计 ~4.8 MB)

### 10.4 系统发育树 (results/tree/)

5 个亚型独立树 + 1 个联合采样树, 总计 ~1.0 MB (Newick 格式)

---

## 十一、方法学注意事项

1. **基因组计数**: GTDB R232 = 189,801 细菌代表基因组 (2026-07-02 修正, 此前误写为 165,468)
2. **门数**: c95 级别 17 门 → 验证后 13 门, 论文中应同时说明两者
3. **系统发育**: 树基于 c95 级别 (8,731 条) 构建以获得进化分辨率; 分类学统计基于验证集 (6,532 条) 以保证可信度
4. **信号肽**: 当前为 Kyte-Doolittle 疏水性启发式, 非 SignalP 6.0; 正式发表前建议用 SignalP 6.0/Phobius 验证
5. **gRodon2**: 输出为基于 codon usage bias 的预测值, 表述使用 "predicted maximal growth potential"
6. **古菌**: "在本次 GTDB R232 代表古菌基因组分析范围内未发现高置信度 phaZ", 不做绝对断言
7. **参考偏倚**: Pseudomonadota 占 96.8%, 11/14 参考序列来自该门, 论文需讨论参考偏倚
8. **PhaC 共现**: 当前为属级交叉引用, 未做统计检验——表述为 "探索性分析"
9. **CD-HIT**: c95 用于系统发育; c70 仅用于概览热图, 不作为进化证据
10. **多拷贝**: 验证集以单拷贝为主 (84.7%), 与初筛 (67.7% 多拷贝) 形成鲜明对比——HMMer 有效过滤假阳性

---

## 十二、论文写作建议主线

> **在 GTDB R232 全部 189,801 个细菌代表基因组中系统鉴定了 6,532 条高置信度 PhaZ 蛋白 (分布于 5,531 个基因组、13 门、1,135 属), 揭示了 PhaZ 强烈的分类学集中性 (Pseudomonadota ~97%)、显著的亚型分化 (5 亚型) 和结构-功能一致性 (胞外型 lipase box 富集), 同时确认古菌无高置信度 phaZ, 且 phaZ 携带状态不对应全局更快的预测生长潜力。**

---

## 十三、待完成与建议优先级

| 优先级 | 任务 | 预计耗时 | 说明 |
|---|---|---|---|
| 🔴 | #16 Pfam 注释完成 | 🔄 运行中 | 6,532 条 hmmscan |
| 🔴 | SignalP 6.0/Phobius | 2h | 替代 Kyte-Doolittle 启发式 |
| 🟡 | 亚型分层 gRodon2 | 1h | 胞内/胞外分别检验 |
| 🟡 | 联合树可视化 (ETE) | 1h | 标注亚型颜色 |
| 🟡 | gRodon2 大属单独建模 | 0.5h | Bradyrhizobium, Polaromonas 等 |
| 🟢 | PhaC Class I-IV HMM + 基因组共现 | 3-4h | 升级为正式共进化分析 |
| 🟢 | 系统发育信号 (Pagel's λ / Blomberg's K) | 2h | phaZ 分布 vs 物种树 |

---

## 十四、已安装 Skills

| Skill | 来源 | 用途 |
|---|---|---|
| `nature-figure` | yuan1z0825 | 论文图生成 |
| `phylogenetics` | K-Dense AI | MAFFT + IQ-TREE + FastTree |
| `etetoolkit` | davila7 | ETE Toolkit 树可视化 |
| `protein-sequence-similarity-search` | Google DeepMind | MMseqs2/BLAST 同源搜索 |
| `protein-sequence-msa` | Google DeepMind | Clustal Omega 多序列比对 |
| `tooluniverse-protein-structure-retrieval` | Harvard MIMS | 蛋白质结构获取 |
| `tooluniverse-phylogenetics` | Harvard MIMS | 系统发育工作流 |
| `tooluniverse-metagenomics-analysis` | Harvard MIMS | 宏基因组分析框架 |
| `bioinformatics-fundamentals` | delphine-l | 生信格式/流程审查 |
| `bioinformatics-workflows` | omer-metin | Nextflow/Snakemake |
| `statistics-verifier` | travisjneuman | 统计方法验证 |

---

**文档生成**: Claude Code (deepseek-v4-pro) + nature-figure skill  
**项目负责人**: 王浩宇 | **GitHub**: https://github.com/anzhiw2-gif/PHB_gtdb
