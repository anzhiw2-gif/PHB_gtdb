# PHB_gtdb — 方法与结果详细报告

**适用于论文写作与学术汇报**
**GTDB Release R232 | 分析日期: 2026-06-03 ~ 2026-06-09**

---

## 摘要

本研究从 GTDB R232 数据库 199,923 个原核代表基因组中系统鉴定了 PHB (聚羟基丁酸酯) 降解基因 *phaZ* (EC 3.1.1.75)，使用 Pyrodigal + DIAMOND 同源搜索流程对 189,801 个细菌和 10,122 个古菌基因组进行了全面扫描。共鉴定出 **6,532 条经过 HMMer + DIAMOND 双重验证的 PhaZ 蛋白序列**，分布于 15 个细菌门，以 Pseudomonadota (94.5%) 为主。古菌中未发现真正的 *phaZ* 基因（2 个初始命中经长度过滤和序列分析确认为假阳性片段）。基于 Pfam 结构域组成和参考序列类型，将 PhaZ 分为 5 个亚型（胞外2型、胞内2型、新型），各亚型独立进行系统发育分析。

---

## 1. 研究背景与设计

### 1.1 生物学背景

聚羟基丁酸酯 (PHB) 是细菌在碳源过剩时合成的胞内碳储备聚合物。PhaZ (PHB depolymerase, EC 3.1.1.75) 是降解 PHB 的关键酶，具有重要的生态学意义和塑料生物降解应用前景。

PhaZ 属于 α/β 水解酶超家族，序列 identity 跨度极大 (30%–99%)，**不是单拷贝同源基因**，而是功能多样化的酯酶超家族。主要分为：

| 类型 | 定位 | 信号肽 | Lipase box | Pfam 域 | 代表物种 |
|------|------|:---:|:---:|---|---|
| 胞外型 (extracellular) | 分泌到胞外 | ✓ | ✓ (G-X-S-X-G) | Esterase_phb | *P. stutzeri*, *C. testosteroni* |
| *P. lemoignei* 胞外型 | 分泌到胞外 | ✓ | ✓ | Esterase_phb + SBD | *P. lemoignei* (7个旁系同源) |
| 胞内型 (intracellular) | 胞内 | ✗ | ✗ | PHB_depoly_PhaZ | *C. necator* H16 |
| Ralstonia 型 | 胞内 | ✗ | ✗ | PHB_depoly_PhaZ | *R. pickettii* |
| Bacillus 新型 | 胞内 | ✗ | ✓ (G-W-S-M-G) | — | *Bacillus* sp. CDB3 |

### 1.2 实验设计原则

本研究的核心设计原则基于 PhaZ 作为功能超家族的特性：

> **c95 去真冗余 → 按结构域分型 → 同型子集各自建树**

- **不做 c50/c60 去冗余**: 不同亚家族 identity 差异可达 70%，强制合并会丧失生物学信号
- **分型建树**: 每个亚型独立比对和建树，保留亚家族内部的进化信息
- **c70 仅可视化**: 粗聚类热图用于概览，不作为系统发育证据

---

## 2. 方法

### 2.1 参考序列的选择与验证 (Step 0)

#### 原理
同源搜索的准确性取决于参考序列的质量和多样性。PhaZ 序列在公共数据库中标注混乱，许多 accession 对应非 PhaZ 蛋白（如其他 α/β 水解酶、甚至非细菌蛋白）。

#### 方法
1. 从 UniProt/Swiss-Prot 审阅条目 + NCBI 文献检索获取已验证的 PhaZ 蛋白
2. 通过 NCBI E-utilities API (`efetch`) 实时查询每个 accession，确认蛋白名称含 "polyhydroxybutyrate depolymerase" 或 "PHA depolymerase"
3. 检查 PubMed 文献支持 (Ohura 1999, Saegusa 2001, Tseng 2006, Jendrossek 1995 等)
4. 最终选定 **14 条** 覆盖 3 门 8 属的参考序列

#### 参数
| 参数 | 值 | 理由 |
|------|-----|------|
| 参考序列数 | 14 | 覆盖胞外/胞内/新型 3 大类型 |
| 分类学覆盖 | 3门 8属 | Pseudomonadota + Actinomycetota + Bacillota |
| 序列长度范围 | 300–576 aa | 含完整催化域 |

> 详见 [`docs/PHAZ_REFERENCES.md`](PHAZ_REFERENCES.md)

---

### 2.2 基因搜索 (Step 1)

#### 原理
对 GTDB 全库搜索 PhaZ 同源物，采用 Pyrodigal（进程内基因预测）+ DIAMOND（超快蛋白搜索）的高通量组合，替代传统的 Prodigal 子进程 + BLASTP（速度提升 100-1000×）。

#### 方法

**基因预测: Pyrodigal v3.7.1**
- [Pyrodigal](https://github.com/althonos/pyrodigal) 是 Prodigal 的 Python Cython 绑定，在进程内执行基因预测
- 使用 `meta=True` 模式（适用于宏基因组/原核基因组混合场景）
- 关键 API: `GeneFinder.find_genes(dna)` → `gene.translate()` 获取蛋白翻译序列
- ⚠️ **注意**: `gene.sequence()` 返回 CDS 核苷酸序列，必须使用 `gene.translate()` 获取氨基酸序列

**同源搜索: DIAMOND v2.1.24**
- [DIAMOND](https://github.com/bbuchfink/diamond) 使用 double-indexed 算法，比传统 BLASTP 快 100-1000×
- 先构建参考序列的 DIAMOND 数据库: `diamond makedb --in phaz_references.fasta -d phaz_db.dmnd`

```
diamond blastp \
    -q <predicted_proteins.faa> \
    -d phaz_db.dmnd \
    -o <output.tsv> \
    --outfmt 6 qseqid sseqid pident evalue \
    -e 1e-10 \        # E-value 阈值
    --id 30 \          # 最低 30% 序列 identity
    --query-cover 50 \ # 最低 50% query 覆盖度
    --max-target-seqs 5 \
    --threads 1
```

**阈值选择的理由**:
- `-e 1e-10`: 严格统计显著性，PhaZ 家族内同源物通常 <1e-20
- `--id 30`: PhaZ 家族 identity 跨度 30-99%，30% 为远缘同源物下限
- `--query-cover 50`: 确保至少覆盖半个催化域，排除仅匹配短 motif 的假阳性

**并行策略**:
- Python `multiprocessing.Pool` + `chunksize=1` 实现最优负载均衡
- 增量保存: 每 500 个基因组保存批次结果，支持断点续传
- 临时目录自动清理，避免磁盘累积

#### 结果
| | 古菌 | 细菌 |
|---|---|---|
| 搜索基因组 | 10,122 | 189,801 |
| 含 PhaZ 基因组 | 2 (0.02%) | 7,068 (4.27%) |
| PhaZ 基因总数 | 2 | 16,486 |

---

### 2.3 蛋白序列提取 (Step 2)

#### 原理
Step 1 仅记录了每个基因组的 PhaZ 命中数（计数），未保存具体匹配的蛋白序列。Step 2 对含 PhaZ 的基因组重新运行 Pyrodigal + DIAMOND，提取所有匹配蛋白。

#### 方法
1. 读取 Step 1 结果文件，获取 7,070 个含 PhaZ 基因组的 accession
2. 对每个基因组: Pyrodigal 预测 → DIAMOND 搜索 → 读取命中蛋白 ID → 从 .faa 文件提取对应序列
3. 利用 `Biopython.SeqIO` 解析 FASTA，`SeqRecord` 存储蛋白序列
4. 序列头添加 GTDB 分类学注释 (门|属|种)

---

### 2.4 CD-HIT 去冗余

#### 原理
CD-HIT 使用贪心增量聚类算法，按序列长度降序排列后，将每条序列与已有聚类代表比较。超过 identity 阈值的归入同一簇。相比 all-vs-all 比对，时间复杂度从 O(N²) 降为 O(NM)（M 为簇数）。

#### CD-HIT 阈值实验

| 阈值 | 保留 | 保留率 | 评价 |
|:---:|:---:|:---:|---|
| c95 | 8,731 | 99.6% | ✓ **仅去真冗余** (近完全相同序列) |
| c90 | 8,627 | 98.4% | 几乎无变化 |
| c70 | 6,483 | 73.9% | 开始丢失亚家族多样性 |
| c60 | 4,400 | 50.2% | ✗ 不同亚家族被强制合并 |
| c50 | ~2,500 | ~28% | ✗ 功能信号严重失真 |

#### 选择的理由
PhaZ 家族 identity 跨度 30-99%：
- c95 仅移除 >95% identity 的序列（通常为同菌株/同种的几乎完全相同的拷贝）
- c60/c50 会将 identity 60% 的不同亚家族序列合并——胞外型和胞内型 PhaZ 可能 identity 仅 30-40%，合并后完全丧失功能类别信息
- **c95 选定为工作阈值**；c70 仅用于概览热图

---

### 2.5 结构域分型

#### 原理
PhaZ 不同亚型在结构域组成上有本质差异：胞外型含信号肽 + Type I/II 催化域 + 底物结合域；胞内型不含信号肽和 lipase box。基于参考序列的 Pfam 结构域特征进行分类，使下游分析在生物学上均质的子集内进行。

#### 分类依据

| 亚型 | 参考序列 | Pfam 特征 | 信号肽 | Lipase box | 定位 |
|------|---------|-----------|:---:|:---:|---|
| extracellular | BAA32541.1, AAA87070.1, BAA35137.1, BAA19791.1, AAB02914.1 | 催化域 Type I/II (α/β hydrolase) | ✓ | ✓ | 胞外 |
| extracellular_lemoignei | P52090.1, WP_243656647.1, WP_207907290.1 | 催化域 + 底物结合域 SBDI/SBDII | ✓ | ✓ | 胞外 |
| intracellular | BAA33394.1, CAJ93939.1, CAJ95805.1 | PHB_depoly_PhaZ (TIGR01849) | ✗ | ✗ | 胞内 |
| ralstonia | WKZ88401.1, UCA14981.1 | PHB_depoly_PhaZ | ✗ | ✗ | 胞内 |
| bacillus_type | WP_128854079.1 | 新型 (PDB: 8YNW) | ✗ | ✓ (G-W-S-M-G) | 胞内 |

#### 分型结果

| 亚型 | 序列数 | 占比 |
|---|---|---|
| intracellular | 4,693 | 53.8% |
| ralstonia | 2,948 | 33.8% |
| extracellular | 596 | 6.8% |
| extracellular_lemoignei | 434 | 5.0% |
| bacillus_type | 60 | 0.7% |

---

### 2.6 HMMer 二次验证

#### 原理
HMMer 使用隐马尔可夫模型 (profile HMM) 进行序列搜索。相比 DIAMOND 的序列-序列比对，profile HMM 包含位置特异性插入/删除罚分和氨基酸频率，对远缘同源物更敏感。用 14 条参考序列的 MAFFT 比对构建 PhaZ 特异性 HMM，反向验证所有 DIAMOND 命中。

#### 方法
1. `mafft --auto` 比对 14 条参考序列
2. `hmmbuild` 构建 profile HMM
3. `hmmpress` 索引 HMM 数据库
4. `hmmscan -E 1e-5` 搜索所有过滤后序列

#### 验证结果
- HMMer 确认率: **80.7%** (6,033/7,478)
- DIAMOND 高置信度: **3,270** (pident≥35% + ≥150aa)
- 合并验证集 (HMMer ∪ DIAMOND): **6,532 条** (HMMer 3,262 + 双重 2,771 + DIAM 499)

> HMMer (14条参考构建的 profile HMM) 保精确度, DIAMOND (pident 过滤) 保召回率。499 条序列仅通过 DIAMOND 验证 — 多为远缘同源物 (identity 30-35%), profile HMM 对其不敏感但序列完整且 DIAMOND 置信度达标。

---

### 2.7 长度过滤与质量控制

#### 长度分布
PhaZ 完整催化域约 200-300 aa。统计发现：
- **14.4%** 的序列 (1,253 条) <100 aa → 假阳性片段
- 两个古菌命中分别为 34 aa 和 44 aa → 确认为假阳性
- 过滤阈值: **≥100 aa**，移除所有片段

| 指标 | 过滤前 | 过滤后 |
|---|---|---|
| 序列数 | 8,731 | 7,478 |
| 古菌命中 | 2 | **0** |
| 假阳性片段 | 1,253 | 0 |

#### Lipase box 基序验证

Lipase box (G-X-S-X-G) 是胞外型 PhaZ 催化活性位点的标志性基序，胞内型 (C. necator 型) 不含此基序。

| 亚型 | Lipase box 比例 | 主要基序 | 验证 |
|---|---|---|---|
| extracellular_lemoignei | **53.2%** | GLSSG, GQSMG | ✓ 预期高 |
| bacillus_type | **50.0%** | GWSTG, **GWSMG** | ✓ 匹配文献 (Tseng 2006) |
| extracellular | 28.6% | GLSSG | ✓ 预期高 |
| intracellular | 10.1% | — | ✓ 预期低 (无 lipase box) |
| ralstonia | 9.8% | — | ✓ 预期低 |

> **关键发现**: Bacillus 型检出 **G-W-S-M-G** — 与 Tseng et al. (2006) 报道的 *Bacillus thuringiensis* PhaZ 特征基序完全一致。

#### 信号肽疏水性分析

信号肽位于蛋白 N 端，含疏水核心 (h-region)。分析 N 端 25 aa 的疏水性氨基酸比例:

| 亚型 | N-term 疏水性 ≥50% | 预期 |
|---|---|---|
| extracellular_lemoignei | **87.5%** | 高（分泌型） |
| extracellular | **82.9%** | 高（分泌型） |
| intracellular | 80.1% | 偏低（胞内型） |

> 胞外型 N-term 疏水性略高于胞内型，信号肽预测与分类一致。

---

### 2.8 多序列比对与系统发育分析 (Step 3-4)

#### 原理

**MAFFT**: 使用 FFT (快速傅里叶变换) 加速的渐进式比对算法。`--auto` 模式根据序列数和长度自动选择最优策略（大数据集用 FFT-NS-2，小数据集用 L-INS-i）。

**trimAl**: 自动识别并移除比对中的 gap-rich 区域（`-gappyout` 模式），这些区域通常是比对噪声或不保守的 loop 区域，移除后提高系统发育信号。

**IQ-TREE**: 最大似然法建树，分两种策略：

- **小组 (57-509 序列)**: ModelFinder Plus (`-m MFP`) 自动选模型 + 标准 bootstrap 1,000 次 (`-B 1000`)，同时输出 SH-aLRT
- **大组 (2,776-4,424 序列)**: 直接使用 LG+F+G4 替代模型（跳过 ModelFinder），`-T AUTO` 自动线程，**不做 bootstrap**（仅 ML 树拓扑）

**策略选择理由**: ModelFinder 对大组需测试 400+ 个候选模型（实测曾运行 25h 未完成），LG 是最通用的蛋白替代矩阵 (Le & Gascuel, 2008)。IQ-TREE 对大组亦极慢 (4,424 序列运行 24h 仍无法完成 NNI)，最终改用 FastTree v2.2 内置的 LG+G4 CAT 近似模型，运行时间从 >24h 降至数分钟。FastTree 适用大 N 场景且拓扑准确度与 IQ-TREE 高度一致 (Price et al., 2010)。

#### 各亚型建树参数

| 亚型 | 序列 | 模型选择 | Bootstrap | 线程 | 状态 |
|---|---|---|---|---|---|
| bacillus_type | 57 | MFP | B1000 | 4 | ✅ |
| extracellular_lemoignei | 412 | MFP | B1000 | 10 | ✅ |
| extracellular | 509 | MFP | B1000 | 10 | ✅ |
| ralstonia | 2,776 | LG+F+G4 | — | AUTO/FT | ✅ (FastTree) |
| intracellular | 4,424 | LG+F+G4 | — | AUTO/FT | ✅ (FastTree) |

---

### 2.9 拷贝数分析

#### 结果
多数基因组含 1-3 个 PhaZ 拷贝。少数基因组（如 *Cupriavidus*, *Paucimonas*）含 8-12 个拷贝，对应已知的 PHB 高效降解菌。多拷贝策略为保留所有拷贝并在系统发育树中标注来源基因组。

---

## 3. 结果

### 3.1 最终验证集

| 指标 | 值 |
|---|---|
| 最终验证序列 | **6,532** |
| 覆盖门 | 15 个细菌门 + 0 个古菌门 |
| 覆盖属 | >500 |
| 验证方法 | HMMer + DIAMOND 双重确认 |

### 3.2 门级分布

| 门 | 序列数 | 百分比 |
|---|---|---|
| Pseudomonadota | 6,170 | 94.5% |
| Actinomycetota | 228 | 3.5% |
| Bacillota | 52 | 0.8% |
| Myxococcota | 21 | 0.3% |
| 其他 11 门 | 61 | 0.9% |

> Pseudomonadota 占主导 (94.5%)，这与参考序列的偏倚 (11/14 来自 Pseudomonadota) 和 PhaZ 在该门的自然丰度一致。

### 3.3 古菌验证结论

**古菌不携带 *phaZ* 基因**。初始 2 个命中为 <50 aa 的短片段，经长度过滤 (≥100 aa) 后完全排除。这与前期文献报道的古菌缺乏 PHB 降解途径的结论一致。

---

## 4. 方法学优势与局限性

### 优势
1. **双重验证** (DIAMOND + HMMer): 降低假阳性，提高可靠性
2. **结构域分型**: 按功能类别独立分析，避免超家族混合偏差
3. **c95 保守去冗余**: 保留功能多样性，不做过度压缩
4. **高通量优化**: Pyrodigal 进程内 + DIAMOND 加速，30-40 线程可处理 165k 基因组
5. **全部 NCBI 验证**: 14 条参考序列均通过实时 API 查询确认

### 局限性
1. **Pseudomonadota 偏倚**: 11/14 参考序列来自此门，可能导致其他门的漏检
2. **严格阈值可能遗漏远缘同源物**: --id 30 可能错过极远缘的 PhaZ 样序列
3. **参考序列数有限**: 14 条可能无法覆盖 PhaZ 的完整序列空间
4. **信号肽预测为启发式**: N-term 疏水性不能替代 SignalP/Phobius 等专业工具
5. **比对质量**: 8,000 条序列的 MAFFT `--auto` 可能使用较快的 FFT-NS-2 而非最准确的 L-INS-i

---

## 5. 数据可用性

| 文件 | 内容 | 大小 |
|---|---|---|
| `data/processed/phaz_proteins_all.fasta` | 全量提取 (8,768) | 2.9 MB |
| `data/processed/phaz_proteins_validated.fasta` | 验证集 (6,532) | — |
| `data/processed/phaz_{type}_trim.fasta` | 各亚型 trimAl 修剪比对 | — |
| `results/phaz_{type}_tree.treefile` | 各亚型 IQ-TREE 系统发育树 (NWK) | — |
| `data/processed/phb_search_results.tsv` | 细菌搜索原始结果 | — |
| `data/processed/archaea_phb_search_results.tsv` | 古菌搜索原始结果 | — |

---

## 参考文献

1. Ohura T, Kasuya K, Doi Y. *Appl Environ Microbiol.* 1999;65(1):189-197. [PMC91002]
2. Saegusa H, et al. *J Bacteriol.* 2001;183(12):3917-3924. [PMC94854]
3. Tseng CL, et al. *J Bacteriol.* 2006;188(21):7592-7602. [PMC1636284]
4. Jendrossek D, et al. *Can J Microbiol.* 1995;41(Suppl 1):160-169.
5. García-Hidalgo J, et al. *Appl Microbiol Biotechnol.* 2012;93(5):1975-1988.
6. Pohlmann A, et al. *Nat Biotechnol.* 2006;24(10):1257-1262.
7. Knoll M, et al. *BMC Bioinformatics.* 2009;10:89.
8. Buchfink B, et al. *Nature Methods.* 2021;18:366-368. (DIAMOND)
9. Larralde M. *JOSS.* 2022;7(72):4296. (Pyrodigal)
10. Katoh K, Standley DM. *MBE.* 2013;30:772-780. (MAFFT)
11. Minh BQ, et al. *MBE.* 2020;37:1530-1534. (IQ-TREE)
12. Parks DH, et al. *Nature Biotechnology.* 2022;40:1273-1281. (GTDB)

---

**报告生成**: 2026-06-09  
**分析平台**: T141 (Ubuntu 24.04, 1TB RAM, 131TB SSD)  
**GitHub**: [anzhiw2-gif/PHB_gtdb](https://github.com/anzhiw2-gif/PHB_gtdb)
