# PHB_gtdb — 最终分析报告

**日期**: 2026-06-03 ~ 2026-06-10  
**GTDB**: Release R232 (199,923 代表基因组)  
**方法**: Pyrodigal + DIAMOND + HMMer + CD-HIT + MAFFT + IQ-TREE/FastTree

---

## 1. 执行摘要

从 GTDB R232 中系统鉴定了 **6,532 条经验证的 PhaZ 蛋白序列**，分布于 **15 个细菌门**。古菌基因组中 **未发现真正的 phaZ 基因**（2 个初始 DIAMOND 命中确认为 <50 aa 假阳性片段）。PhaZ 按 Pfam 结构域特征分为 5 个亚型，各亚型独立建树。

---

## 2. 数据统计

| 指标 | 值 |
|---|---|
| 搜索细菌基因组 | 189,801 |
| 搜索古菌基因组 | 10,122 |
| DIAMOND 命中 (细菌) | 7,068 (4.27%) |
| DIAMOND 命中 (古菌) | 2 (0.02%) |
| PhaZ 基因总数 | 16,488 |
| 长度过滤 (≥100aa) | 7,478 (移除 1,253 片段) |
| HMMer 验证通过 | 6,033/7,478 (80.7%) |
| DIAMOND 高置信度 | 3,270 (pident≥35% + ≥150aa) |
| **最终验证集** | **6,532** (HMMer + DIAMOND 并集) |
| 覆盖门 | 15 细菌 + 0 古菌 |
| 覆盖属 | 1,135 |

> **最终验证集 = HMMer 验证 ∪ DIAMOND 高置信度**: 6,033 条通过 profile HMM (E≤1e-5, 精确但保守), 499 条额外通过 DIAMOND 高置信度 (pident≥35% + len≥150aa)。两种方法互补, 确保精确度与召回率的平衡。

---

## 3. 古菌结论

**古菌不携带 phaZ 基因**。

2 个初始 DIAMOND 命中: Halovenus sp. (34 aa) 和 UBA73 sp. (44 aa)。PhaZ 完整催化域 >200 aa，这 2 个序列长度远低于阈值。经 HMMer 验证 + ≥100 aa 过滤后，古菌命中数为 **0**。

---

## 4. 亚型分布

| 亚型 | 序列数 | Identity (mean) | 长度 (mean) | Lipase box | 信号 |
|---|---|---|---|---|---|
| intracellular | 3,668 (56.2%) | 37.1% | 204 aa | 10.1% | 胞内 |
| ralstonia | 2,221 (34.0%) | 37.0% | 205 aa | 9.8% | 胞内 |
| extracellular | 404 (6.2%) | 38.0% | 247 aa | **28.6%** | 胞外 |
| extracellular_lemoignei | 234 (3.6%) | 34.9% | 220 aa | **53.2%** | 胞外 |
| bacillus_type | 5 (0.08%) | 33.5% | 170 aa | **50.0%** | 胞内 |

### Lipase box 验证

胞外型 lipase box (G-X-S-X-G) 比例为 28.6-53.2%，显著高于胞内型 (9.8-10.1%)。Bacillus 型检出 **G-W-S-M-G** 基序，与 Tseng et al. (2006) 报道完全一致。

---

## 5. 门级分布

| 门 | 序列数 | 占比 |
|---|---|---|
| Pseudomonadota | 6,323 | 96.8% |
| Actinomycetota | 154 | 2.4% |
| Myxococcota | 14 | 0.2% |
| 其他 10 门 | 41 | 0.6% |

> Pseudomonadota 占主导 (97%)，部分反映参考序列偏倚 (11/14 来自此门)。

---

## 6. 富集属 Top 10

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

---

## 7. 系统发育树

| 组 | 序列 | 工具 | 模型 | 文件 |
|---|---|---|---|---|
| bacillus_type | 57 | IQ-TREE | MFP + B1000 | 6.3K |
| extracellular_lemoignei | 412 | IQ-TREE | MFP + B1000 | 48K |
| extracellular | 509 | IQ-TREE | MFP + B1000 | 58K |
| ralstonia | 2,776 | IQ-TREE | LG+F+G4 | 317K |
| ralstonia | 2,860 | FastTree | LG+G4 | 325K |
| intracellular | 4,424 | FastTree | LG+G4 | 528K |

> 大组 (>2,000 序列) ModelFinder prohibitive (测试 437 个模型 25h 未完成)。改用 FastTree v2.2 LG+G4 近似极大似然，运行时间从 >24h 降至数分钟。IQ-TREE NNI 对 4,424 序列运行 >24h 仍无法完成。

---

## 8. 管线性能

| 步骤 | 处理量 | 耗时 | 资源 |
|---|---|---|---|
| Step 1 细菌搜索 | 189,801 基因组 | 40.5h | 30 线程 |
| Step 1 古菌搜索 | 10,122 基因组 | 2.2h | 10 线程 |
| Step 2 序列提取 | 7,070 基因组 | 13h | 40 线程 |
| CD-HIT c95 | 8,768 | 2min | 30 线程 |
| 长度过滤 + QС | 8,731 → 6,532 | <1min | — |
| HMMer 验证 | 7,478 | 5min | 30 线程 |
| MAFFT ×5 | ~8,000 | ~10min | 8-20 线程 |
| IQ-TREE ×3 (小组) | 57-509 | ~3h | 4-10 线程 |
| FastTree ×2 (大组) | 2,776-4,424 | ~5min | 1 线程 |
| **总计** | | **~59h** | |

---

## 9. 输出文件

```
results/
├── phaz_bacillus_type_tree.treefile         (6.3K)
├── phaz_extracellular_lemoignei_tree.treefile (48K)
├── phaz_extracellular_tree.treefile          (58K)
├── phaz_ralstonia_tree.treefile             (317K)
├── phaz_ralstonia_tree_ft.treefile          (325K)
├── phaz_intracellular_tree_ft.treefile      (528K)

data/processed/
├── phaz_proteins_all.fasta      (8,768 条, ~3MB)
├── phaz_proteins_validated.fasta (6,532 条)
├── phaz_{type}.fasta            (5 亚型)
├── phaz_{type}_trim.fasta       (trimAl 后)
├── phb_search_results.tsv       (7,068 行)
└── archaea_phb_search_results.tsv (2 行, 0 真阳性)

docs/
├── METHODS.md           论文方法
├── RESULTS.md           统计结果
├── PIPELINE.md          管线流程
├── PHAZ_REFERENCES.md   参考序列文档
├── REVIEW.md            评审与改进计划
└── FINAL_REPORT.md      本文件
```

---

## 10. 参考文献

1. Ohura et al. (1999) *AEM* 65:189-197
2. Saegusa et al. (2001) *J Bacteriol* 183:3917-3924
3. Tseng et al. (2006) *J Bacteriol* 188:7592-7602
4. Jendrossek et al. (1995) *Can J Microbiol* 41:160-169
5. García-Hidalgo et al. (2012) *AMB* 93:1975-1988
6. Buchfink et al. (2021) *Nat Methods* 18:366-368
7. Larralde (2022) *JOSS* 7:4296
8. Katoh & Standley (2013) *MBE* 30:772-780
9. Minh et al. (2020) *MBE* 37:1530-1534
10. Price et al. (2010) *PLoS ONE* 5:e9490
11. Parks et al. (2022) *Nat Biotechnol* 40:1273-1281

---

**报告生成**: 2026-06-10 | **平台**: Ubuntu 24.04, 1TB RAM, 131TB SSD | **GitHub**: [anzhiw2-gif/PHB_gtdb](https://github.com/anzhiw2-gif/PHB_gtdb)
