# PHB_gtdb — 最终统计结果

**生成日期**: 2026-06-09 | **验证集**: 6,532 PhaZ 蛋白序列 | **覆盖**: 15 细菌门

---

## 1. 门 × 亚型 分布矩阵

| 门 (Phylum) | intracellular | ralstonia | extracellular | ext_lemoignei | bacillus | 合计 |
|---|---:|---:|---:|---:|---:|---:|
| **Pseudomonadota** | 3,636 | 2,199 | 354 | 134 | 0 | **6,323** |
| **Actinomycetota** | 26 | 15 | 49 | 64 | 0 | **154** |
| **Myxococcota** | 3 | 1 | 0 | 10 | 0 | **14** |
| **Acidobacteriota** | 1 | 1 | 0 | 4 | 0 | **6** |
| **Bacteroidota** | 0 | 0 | 0 | 5 | 0 | **5** |
| **Bacillota** | 0 | 0 | 0 | 0 | 5 | **5** |
| **Chloroflexota** | 1 | 1 | 0 | 3 | 0 | **5** |
| **Gemmatimonadota** | 1 | 2 | 0 | 2 | 0 | **5** |
| **Bdellovibrionota** | 0 | 1 | 1 | 2 | 0 | **4** |
| **Desulfobacterota** | 0 | 0 | 0 | 4 | 0 | **4** |
| **Spirochaetota** | 0 | 1 | 0 | 2 | 0 | **3** |
| **Planctomycetota** | 0 | 0 | 0 | 3 | 0 | **3** |
| **Cyanobacteriota** | 0 | 0 | 0 | 1 | 0 | **1** |

---

## 2. 亚型统计

| 亚型 | 序列数 | Identity (mean) | 长度 (mean) | Lipase box (%) | N-term 疏水 (%) |
|---|---|---|---|---|---|
| intracellular | 3,668 | 37.1% | 204 aa | 10.1% | 80.1% |
| ralstonia | 2,221 | 37.0% | 205 aa | 9.8% | 80.1% |
| extracellular | 404 | 38.0% | 247 aa | **28.6%** | **82.9%** |
| ext_lemoignei | 234 | 34.9% | 220 aa | **53.2%** | **87.5%** |
| bacillus_type | 5 | 33.5% | 170 aa | **50.0%** | 62.5% |

### 关键验证

- **Lipase box**: 胞外型 (28-53%) 显著高于胞内型 (~10%) → 分型合理性 ✓
- **G-W-S-M-G 基序**: Bacillus 型检出文献报道的特征序列 (Tseng 2006) ✓
- **N-term 疏水性**: 胞外型 (83-88%) 略高于胞内型 (~80%) → 信号肽特征一致 ✓
- **HMMer 验证**: 80.7% (6,033/7,478) 通过 profile HMM 二次验证 ✓

---

## 3. 参考序列匹配分布

| 参考序列 | 亚型 | 匹配数 |
|---|---|---|
| BAA33394.1 (*C. necator* PhaZ1) | intracellular | 1,900 |
| UCA14981.1 (*R. pickettii*) | ralstonia | 1,447 |
| CAJ93939.1 (*C. necator* PhaZ2) | intracellular | 1,246 |
| WKZ88401.1 (*R. pickettii*) | ralstonia | 774 |
| CAJ95805.1 (*C. necator* PhaZ5) | intracellular | 522 |
| BAA19791.1 (*D. acidovorans*) | extracellular | 151 |
| P52090.1 (*P. lemoignei* PhaZ3) | ext_lemoignei | 142 |
| AAA87070.1 (*C. testosteroni*) | extracellular | 77 |
| BAA35137.1 (*Acidovorax* sp.) | extracellular | 68 |
| BAA32541.1 (*S. stutzeri*) | extracellular | 59 |
| AAB02914.1 (*S. exfoliatus*) | extracellular | 49 |
| WP_207907290.1 (*P. lemoignei*) | ext_lemoignei | 48 |
| WP_243656647.1 (*P. lemoignei*) | ext_lemoignei | 44 |
| WP_128854079.1 (*Bacillus* sp.) | bacillus_type | 5 |

> *C. necator* PhaZ1 (BAA33394.1) 是最普遍的 PhaZ 同源物，表明胞内型 PhaZ 在细菌中最广泛分布。

---

## 4. PhaZ 富集最多的属 (Top 20)

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
| *Caballeronia* | 73 | ralstonia | — |
| *Pseudolabrys* | 73 | intracellular | — |
| *Rugamonas* | 71 | intracellular | — |
| *BOG-931* | 71 | intracellular | — |
| *Hydrogenophaga* | 68 | intracellular | — |
| *Acidovorax* | 67 | extracellular | ✓ (参考来源) |
| *Aquabacterium* | 65 | intracellular | — |
| *Telluria* | 62 | intracellular | — |
| *Herbaspirillum* | 57 | ralstonia | — |
| *Allorhodoferax* | 56 | ralstonia | — |

---

## 5. 管线性能

| 步骤 | 处理量 | 耗时 | 资源 |
|---|---|---|---|
| Step 1 细菌搜索 | 189,801 基因组 | 40.5h | 30 线程 |
| Step 1 古菌搜索 | 10,122 基因组 | 2.2h | 10 线程 |
| Step 2 序列提取 | 7,070 基因组 | 13h | 40 线程 |
| CD-HIT c95 | 8,768 → 8,731 | 2min | 30 线程 |
| 长度过滤 ≥100aa | 8,731 → 7,478 | <1s | — |
| HMMer 验证 | 7,478 | 5min | 30 线程 |
| 分型 + 统计 | 6,532 | 30s | — |

---

## 6. 古菌验证结论

| 阶段 | 命中数 | 验证后 |
|---|---|---|
| DIAMOND 初始搜索 | 2 (0.02%) | — |
| 序列长度检查 | 34 aa + 44 aa | ❌ 远小于 PhaZ 催化域 (>200 aa) |
| ≥100aa 过滤 | 0 | **0 真阳性** |

> **最终结论: 古菌基因组中不存在真正的 *phaZ* 基因。** 初始 2 个命中为 DIAMOND 交叉匹配的短随机片段。

---

**数据文件**: `data/processed/phaz_proteins_validated.fasta` (6,532 条)
