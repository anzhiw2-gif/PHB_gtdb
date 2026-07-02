# PHB_gtdb 项目汇报总结

**建议汇报题目**: GTDB R232 中 PHB 降解关键酶 PhaZ 的系统鉴定、分型与系统发育分析  
**项目对象**: PHB/PHA depolymerase，即 *phaZ* 编码的聚羟基丁酸酯降解相关蛋白  
**分析数据**: GTDB Release R232 代表基因组  
**运行平台**: T141 服务器，项目路径 `/home/data/haoyu/PHB_gtdb`  
**核心流程**: Pyrodigal + DIAMOND + CD-HIT + HMMer + MAFFT/trimAl + IQ-TREE/FastTree

---

## 1. 一句话总结

本项目在 GTDB R232 代表基因组中系统搜索 PHB 降解关键酶 PhaZ，最终获得 **6,532 条高置信度 PhaZ 蛋白序列**，覆盖 **15 个细菌门、1,135 个属**；在本次分析范围内，古菌中未发现高置信度 *phaZ* 基因。结果显示 PhaZ 高度富集于 **Pseudomonadota**，并以 **胞内型 PhaZ** 为主，提示 PHB 代谢相关酶在细菌中具有明显的分类学集中性和亚型分化。

---

## 2. 研究背景与科学问题

PHB，聚羟基丁酸酯，是许多微生物合成和储存的一类可降解聚合物，也是 PHA 生物塑料的重要类型。PhaZ 是 PHB/PHA 降解过程中的关键酶，决定微生物是否具备利用或降解 PHB 的潜力。

本项目主要回答以下问题：

1. 在 GTDB R232 代表基因组中，哪些微生物携带 PhaZ？
2. PhaZ 是否主要集中在特定细菌类群中？
3. 古菌是否存在可信的 *phaZ* 基因？
4. PhaZ 能否根据参考序列、结构特征和保守基序分成稳定亚型？
5. 不同 PhaZ 亚型是否应该分开进行系统发育分析？

---

## 3. 数据来源与参考序列

### 3.1 基因组数据

| 数据类型 | 数量 | 说明 |
|---|---:|---|
| GTDB R232 代表基因组总量 | 199,923 | 项目使用的数据库背景 |
| 本次搜索的细菌基因组 | 189,801 | 主要搜索对象 |
| 本次搜索的古菌基因组 | 10,122 | 用于验证古菌是否携带 *phaZ* |

### 3.2 PhaZ 参考序列

项目使用 **14 条经 NCBI/UniProt 与文献验证的 PhaZ 参考序列**，覆盖胞外分泌型、胞内型、Ralstonia 相关型和 Bacillus 新型 PhaZ。

| 参考类型 | 代表来源 | 生物学特征 |
|---|---|---|
| extracellular | *Stutzerimonas*, *Comamonas*, *Acidovorax*, *Delftia* 等 | 多为胞外分泌型，常见 G-X-S-X-G lipase box |
| extracellular_lemoignei | *Paucimonas lemoignei* 多个旁系同源 PhaZ | 胞外型，具有较明显的分支特征 |
| intracellular | *Cupriavidus necator* H16 PhaZ1/2/5 | 胞内 PHB 颗粒降解相关 |
| ralstonia | *Ralstonia pickettii* 相关序列 | 胞内型变体 |
| bacillus_type | *Bacillus* sp. CDB3 | 新型 PhaZ，具有 G-W-S-M-G 特征基序 |

---

## 4. 分析流程

项目的整体分析逻辑是：先在大规模基因组中搜索候选 PhaZ，再通过去冗余、长度过滤、HMMer 验证和 DIAMOND 高置信度标准得到最终集合，最后进行亚型划分、系统发育分析和可视化。

```text
GTDB R232 代表基因组
  -> Pyrodigal 预测蛋白编码基因
  -> DIAMOND blastp 搜索 PhaZ 同源序列
  -> 提取候选 PhaZ 蛋白
  -> CD-HIT c95 去除近完全重复序列
  -> 长度过滤，去除短片段
  -> HMMer + DIAMOND 高置信度联合验证
  -> 根据参考序列与结构特征划分 5 个亚型
  -> 各亚型分别进行 MAFFT 比对和 trimAl 修剪
  -> IQ-TREE 或 FastTree 构建系统发育树
  -> 统计分布、保守基序和可视化结果
```

### 为什么按亚型分别建树？

PhaZ 不是单一高度保守蛋白，而是属于功能多样的 alpha/beta hydrolase 相关超家族，不同亚型之间序列相似度和结构域组成差异较大。胞外型、胞内型、Ralstonia 型和 Bacillus 新型 PhaZ 在定位、保守基序和长度上都不同。如果把全部序列混合建一棵总树，容易造成不同功能类别之间的错误比较。

因此本项目采用：

```text
c95 去冗余 -> 亚型分型 -> 各亚型独立比对和建树
```

这一策略更符合 PhaZ 的生物学差异，也能让系统发育结果更容易解释。

---

## 5. 主要结果总览

### 5.1 从初始搜索到最终验证集

| 阶段 | 数量 | 说明 |
|---|---:|---|
| 细菌 DIAMOND 命中基因组 | 7,068 | 占搜索细菌基因组约 4.27% |
| 细菌 PhaZ 初始命中数 | 16,486 | Step 1 初始同源搜索结果 |
| 古菌 DIAMOND 初始命中 | 2 | 后续判定为短片段假阳性 |
| 提取的候选 PhaZ 蛋白 | 8,768 | 全量候选蛋白序列 |
| CD-HIT c95 后 | 8,731 | 去除近完全重复序列 |
| 长度过滤后，>=100 aa | 7,478 | 去除 1,253 条短片段 |
| HMMer 验证通过 | 6,033 | profile HMM，E <= 1e-5 |
| DIAMOND 高置信度补充 | 499 | pident >=35% 且长度 >=150 aa |
| 最终验证集 | **6,532** | HMMer 与 DIAMOND 高置信度结果的并集 |

### 5.2 最终结果一句话概括

最终得到的 **6,532 条 PhaZ 序列** 不是单纯 DIAMOND 命中结果，而是经过 **去冗余、长度过滤、HMMer 验证和高置信度 DIAMOND 补充** 后得到的验证集合，因此更适合作为后续分布统计和系统发育分析的基础。

---

## 6. 核心发现一：古菌中未发现高置信度 *phaZ*

古菌搜索只得到 2 个 DIAMOND 初始命中，但对应序列长度仅约 **34 aa** 和 **44 aa**，远低于完整 PhaZ 催化域通常需要的长度。经过长度过滤和 HMMer 验证后，古菌中没有序列进入最终验证集。

| 阶段 | 古菌命中数 | 解释 |
|---|---:|---|
| DIAMOND 初始命中 | 2 | 仅为初始同源搜索结果 |
| 长度过滤后 | 0 | 两条序列均为短片段 |
| 最终验证集 | 0 | 未发现高置信度古菌 PhaZ |

**汇报表述建议**: 在本次 GTDB R232 代表古菌基因组分析范围内，未发现可信的 *phaZ* 基因；初始 2 个命中更可能是短片段交叉匹配或随机同源片段，而不是完整功能基因。

---

## 7. 核心发现二：PhaZ 高度集中于 Pseudomonadota

最终验证集覆盖 15 个细菌门，但分布极不均匀，其中 **Pseudomonadota 占 96.8%**，是最主要的 PhaZ 来源类群。

| 门 | PhaZ 序列数 | 占比 |
|---|---:|---:|
| Pseudomonadota | 6,323 | 96.8% |
| Actinomycetota | 154 | 2.4% |
| Myxococcota | 14 | 0.2% |
| 其他门合计 | 41 | 0.6% |

这一结果说明 PhaZ 在 GTDB 代表细菌基因组中具有明显的分类学集中性。需要注意的是，当前参考序列也主要来自 Pseudomonadota，因此该结果既反映真实富集趋势，也可能受到参考序列来源偏倚的影响。

---

## 8. 核心发现三：胞内型 PhaZ 是主导类型

根据参考序列和结构特征，最终 PhaZ 验证集分为 5 个亚型。

| 亚型 | 序列数 | 占比 | 主要解释 |
|---|---:|---:|---|
| intracellular | 3,668 | 56.2% | 最主要的胞内型 PhaZ |
| ralstonia | 2,221 | 34.0% | 胞内型相关变体 |
| extracellular | 404 | 6.2% | 胞外分泌型 |
| extracellular_lemoignei | 234 | 3.6% | *P. lemoignei* 相关胞外型 |
| bacillus_type | 5 | 0.08% | Bacillus 新型 PhaZ |

`intracellular` 与 `ralstonia` 两类合计超过 **90%**。这说明在 GTDB 代表基因组中，PhaZ 更主要地表现为胞内 PHB 储存颗粒代谢相关酶，而不是胞外分泌型降解酶。

---

## 9. 核心发现四：Lipase Box 支持当前亚型划分

Lipase box，尤其是 G-X-S-X-G，是许多水解酶和胞外型 PhaZ 的重要催化相关基序。本项目发现胞外型 PhaZ 的 lipase box 检出比例明显高于胞内型。

| 亚型 | Lipase box 检出比例 | 解释 |
|---|---:|---|
| extracellular_lemoignei | 53.2% | 胞外型特征明显 |
| extracellular | 28.6% | 高于胞内型 |
| bacillus_type | 50.0% | 检出 G-W-S-M-G 特征基序 |
| intracellular | 10.1% | 较低，符合胞内型预期 |
| ralstonia | 9.8% | 较低，符合胞内型预期 |

这个结果从保守基序层面支持当前分型策略：胞外型与 Bacillus 新型更容易保留典型催化基序，而胞内型 PhaZ 的序列特征相对不同。

---

## 10. 核心发现五：富集属可作为后续实验筛选对象

PhaZ 最富集的属多属于已知或潜在 PHB 代谢相关类群，说明本项目结果可以为后续菌株筛选提供候选名单。

| 属 | PhaZ 序列数 | 主导亚型 | 汇报解释 |
|---|---:|---|---|
| *Paraburkholderia* | 319 | ralstonia | PhaZ 富集最明显 |
| *Bradyrhizobium* | 239 | intracellular | 胞内型为主 |
| *Cupriavidus* | 148 | intracellular | 经典 PHB 代谢参考来源 |
| *Variovorax* | 136 | extracellular | 胞外型相关 |
| *Polaromonas* | 109 | intracellular | 潜在候选类群 |
| *Rubrivivax* | 107 | intracellular | 已知相关类群 |
| *Polynucleobacter* | 98 | intracellular | 值得进一步验证 |
| *Mesorhizobium* | 96 | intracellular | 已知相关类群 |
| *Stenotrophomonas* | 93 | intracellular | 潜在候选类群 |
| *Burkholderia* | 90 | ralstonia | 已知相关类群 |

---

## 11. 系统发育分析结果

本项目对 5 个 PhaZ 亚型分别进行多序列比对和建树。小规模亚型使用 IQ-TREE 进行模型选择和 bootstrap 分析；大规模亚型由于序列数量达到数千条，使用 FastTree 作为近似极大似然方案。

| 亚型 | trim 后序列数 | 建树策略 |
|---|---:|---|
| bacillus_type | 57 | IQ-TREE + ModelFinder + bootstrap |
| extracellular | 509 | IQ-TREE + ModelFinder + bootstrap |
| extracellular_lemoignei | 412 | IQ-TREE + ModelFinder + bootstrap |
| ralstonia | 2,776 | IQ-TREE / FastTree 辅助 |
| intracellular | 4,424 | FastTree，大规模近似 ML |

**汇报重点**: 不是所有 PhaZ 强行混合建树，而是根据亚型分别建树。这一点是本项目方法设计中比较重要的合理性来源。

---

## 12. 已生成图表与汇报用途

建议汇报时优先使用以下图表：

| 图表文件 | 建议用途 |
|---|---|
| `results/figures/pipeline_timeline.png` | 展示整体流程和运行耗时 |
| `results/figures/phylum_subtype_heatmap.png` | 展示门水平与亚型的联合分布，是核心结果图 |
| `results/figures/subtype_pie.png` | 展示 5 个 PhaZ 亚型的比例 |
| `results/figures/phylum_barplot.png` | 展示 PhaZ 在不同门中的集中分布 |
| `results/figures/top_genera.png` | 展示 PhaZ 富集属 |
| `results/figures/lipase_box_barplot.png` | 展示 lipase box 对分型的支持 |
| `results/figures/reference_matches.png` | 展示参考序列匹配来源 |
| `results/figures/length_distribution.png` | 展示长度过滤和序列质量控制的合理性 |

如果时间有限，建议只放 4 张核心图：

1. `pipeline_timeline.png`
2. `phylum_subtype_heatmap.png`
3. `subtype_pie.png`
4. `top_genera.png`

---

## 13. 推荐 PPT 结构

### 第 1 页：题目与研究目标

说明本项目目标是在 GTDB R232 中系统鉴定 PHB 降解关键酶 PhaZ，并分析其分类分布、功能亚型和系统发育关系。

### 第 2 页：研究背景

介绍 PHB/PHA 生物塑料、PhaZ 的降解作用，以及研究 PhaZ 对塑料生物降解和微生物碳循环的意义。

### 第 3 页：数据与参考序列

展示 GTDB R232、细菌和古菌基因组数量，以及 14 条经验证 PhaZ 参考序列。

### 第 4 页：分析流程

用流程图展示 Pyrodigal、DIAMOND、CD-HIT、HMMer、MAFFT、IQ-TREE/FastTree 的整体管线。

### 第 5 页：数据筛选过程

用漏斗图或表格展示从 16,486 个细菌初始命中到 6,532 条最终验证序列的筛选过程。

### 第 6 页：古菌验证结果

突出古菌只有 2 个短片段初始命中，最终没有高置信度 *phaZ*。

### 第 7 页：门水平分布

展示 Pseudomonadota 占主导，Actinomycetota 次之，其他门数量较少。

### 第 8 页：亚型组成

展示 intracellular 与 ralstonia 两类胞内型 PhaZ 占主导，胞外型比例较低。

### 第 9 页：Lipase box 与分型验证

说明胞外型 lipase box 检出率更高，支持当前亚型划分。

### 第 10 页：富集属与潜在应用

展示 *Paraburkholderia*、*Bradyrhizobium*、*Cupriavidus*、*Variovorax* 等富集属，提出后续菌株筛选方向。

### 第 11 页：系统发育分析

说明为什么不同亚型需要分开建树，以及大规模亚型为什么使用 FastTree。

### 第 12 页：结论与下一步

总结核心发现、项目优势、局限性和后续验证计划。

---

## 14. 汇报结论

可以在最后一页使用以下结论：

1. 本项目在 GTDB R232 中系统鉴定出 **6,532 条高置信度 PhaZ 蛋白序列**，覆盖 **15 个细菌门和 1,135 个属**。
2. 在本次分析范围内，古菌中未发现高置信度 *phaZ* 基因；初始 2 个古菌命中均为短片段假阳性。
3. PhaZ 在细菌中高度集中于 **Pseudomonadota**，尤其富集于 *Paraburkholderia*、*Bradyrhizobium*、*Cupriavidus* 和 *Variovorax* 等属。
4. PhaZ 以胞内型为主，`intracellular` 与 `ralstonia` 两类合计超过 **90%**。
5. Lipase box 分布支持当前亚型划分，胞外型和 Bacillus 新型比胞内型更容易保留典型催化基序。
6. PhaZ 是一个功能多样的酶类群，因此采用 c95 去冗余、亚型分型和分亚型建树，比混合建一棵总树更合理。

---

## 15. 项目亮点

- 使用 GTDB R232 代表基因组，数据覆盖范围较大。
- 同时分析细菌和古菌，能够回答古菌是否携带 *phaZ* 的问题。
- 参考序列经过文献和数据库验证，降低错误注释带来的影响。
- 使用 DIAMOND 与 HMMer 联合验证，兼顾搜索召回率和结果可信度。
- 采用 c95 去除近完全重复序列，同时保留 PhaZ 的亚型多样性。
- 按 PhaZ 亚型分别建树，避免不同功能类型混合造成系统发育解释偏差。
- 结果不仅给出分布规律，也提供了后续实验筛选的候选菌属和候选序列。

---

## 16. 局限性与下一步工作

### 局限性

1. 参考序列主要来自 Pseudomonadota，可能导致搜索结果对该门更敏感。
2. 当前信号肽和胞外定位主要通过序列特征间接判断，仍需 SignalP 或 Phobius 进一步验证。
3. HMMer profile 依赖现有参考序列组成，可能对远缘 PhaZ 不够敏感。
4. 大规模亚型树使用 FastTree，是计算效率和模型精度之间的折中。
5. 当前结果是计算预测结果，关键候选菌株仍需要实验验证其 PHB 降解能力。

### 下一步工作

1. 扩充非 Pseudomonadota 来源的 PhaZ 参考序列，降低参考偏倚。
2. 使用 SignalP/Phobius 验证胞外型 PhaZ 的信号肽特征。
3. 使用 InterProScan/Pfam 系统注释最终验证序列的结构域组成。
4. 对少数非 Pseudomonadota 门中的命中进行人工核查。
5. 选择富集属或代表性候选菌株进行 PHB 降解实验验证。
6. 对胞内型大类进一步分析拷贝数、物种分布和潜在水平基因转移事件。

---

## 17. 口头汇报稿参考

本项目围绕 PHB 降解关键酶 PhaZ，在 GTDB R232 代表基因组中进行了系统搜索、验证、分型和系统发育分析。我们使用 14 条经过文献和数据库验证的 PhaZ 参考序列，结合 Pyrodigal 基因预测和 DIAMOND 同源搜索，对 189,801 个细菌基因组和 10,122 个古菌基因组进行了筛查。经过 CD-HIT c95 去冗余、长度过滤、HMMer 验证以及 DIAMOND 高置信度补充后，最终获得 6,532 条高置信度 PhaZ 蛋白序列。

结果显示，PhaZ 主要分布于细菌中，覆盖 15 个细菌门和 1,135 个属；在古菌中仅发现 2 个很短的初始命中，经过过滤和验证后均未进入最终结果，因此在本次分析范围内没有发现可信的古菌 *phaZ*。从分类分布看，PhaZ 高度集中于 Pseudomonadota，占最终验证集的 96.8%，其中 *Paraburkholderia*、*Bradyrhizobium*、*Cupriavidus* 和 *Variovorax* 等属最为富集。

从功能亚型看，PhaZ 可以分为 intracellular、ralstonia、extracellular、extracellular_lemoignei 和 bacillus_type 五类。其中 intracellular 和 ralstonia 两类胞内型合计超过 90%，说明在 GTDB 代表基因组中，PhaZ 更多与胞内 PHB 储存颗粒代谢相关，而胞外分泌型比例较低。Lipase box 的分布也支持这一分型结果：胞外型和 Bacillus 新型中 lipase box 检出率更高，而胞内型明显较低。

在系统发育分析中，我们没有把所有 PhaZ 序列混合建树，而是根据亚型分别进行比对和建树。这样做可以避免胞外型、胞内型和新型 PhaZ 之间因序列差异过大而造成解释偏差。总体来看，本项目建立了一套适合大规模基因组筛查的 PhaZ 分析流程，获得了一批高置信度候选序列和潜在 PHB 降解菌属，为后续塑料生物降解菌株筛选和功能实验提供了数据基础。

---

## 18. 汇报时需要注意的表述

- 建议说“在本次 GTDB R232 代表基因组分析范围内未发现高置信度古菌 *phaZ*”，不要绝对说“所有古菌都没有 *phaZ*”。
- 建议说“Pseudomonadota 高度富集”，同时承认参考序列主要来自该门，存在一定参考偏倚。
- 建议说“最终验证集为 HMMer 与 DIAMOND 高置信度结果的并集”，避免听众误解为只用了单一搜索标准。
- 建议强调“分亚型建树”的合理性，这是本项目方法设计中的关键点。
- 汇报图表优先展示流程、门水平热图、亚型比例和富集属，不必在主汇报中展示所有中间文件。
