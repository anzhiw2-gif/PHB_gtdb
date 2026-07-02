# 数据账本：GTDB R232 基因组扫描与过滤追踪

**日期**: 2026-07-02 | **版本**: v1.0

---

## 1. GTDB 官方数据

| 指标 | 数值 | 来源 |
|---|---|---|
| GTDB Release | R232 | gtdb.ecogenomic.org |
| Total genome accessions | 878,998 (bac120) + 23,343 (ar53) | `*_taxonomy_r232.tsv` |
| Species clusters | 199,923 | GTDB official |
| Bacterial species clusters | **189,801** | GTDB official |
| Archaeal species clusters | **10,122** | GTDB official |

---

## 2. 本项目扫描的基因组

| 指标 | 数值 | 来源 |
|---|---|---|
| 数据库文件数 | 199,923 `*_genomic.fna.gz` | `gtdb_genomes_reps_r232/database/` |
| 扫描总数 | 199,923 | pipeline log |
| 细菌 (Domain=Bacteria) | **189,801** | domain_map 分类 |
| 古菌 (Domain=Archaea) | **10,122** | domain_map 分类 |
| 未知域 | 0 | 全部成功分类 |

> ⚠️ **之前 README 误写为 165,468**。该数字可能是早期 GTDB 版本的计数（如 R220），或对基因组列表过滤后的子集。已于 2026-07-02 修正为 189,801。

---

## 3. 管道各级过滤

| 步骤 | 序列数 | 基因组数 | 移除项 |
|---|---|---|---|
| DIAMOND 原始命中 | 16,486 hits | 7,068 | — |
| 蛋白序列提取 | 8,767 | 7,068 | 去重叠 CDS |
| CD-HIT c95 | 8,731 | 7,051 | 36 条 (0.4%) |
| 长度 ≥100 aa | 7,478 | 6,224 | 1,253 短片段 (含古菌 2 条) |
| **HMMer + DIAMOND 验证** | **6,532** | **5,531** | 946 条未通过验证 |
| trimAl 修剪 | 8,178 | — | 基于 c95 级别比对 |

---

## 4. 门数变化

| 级别 | 门数 | 新增/移除的门 |
|---|---|---|
| CD-HIT c95 | **17** | 含 Halobacteriota, Nanobdellota, Thermotogota, Verrucomicrobiota |
| 验证集 | **13** | 移除 4 个古菌相关/稀少门 |

> 论文中应明确：**"初始 c95 级别覆盖 17 门，经严格验证后最终高置信度集包含 13 门"**。

---

## 5. c95 vs 验证集

| 亚型 | c95 序列 | c95 基因组 | 验证序列 | 验证基因组 | 验证率 |
|---|---|---|---|---|---|
| intracellular.cupriavidus_like | 4,693 | 4,200 | 3,668 | 3,358 | 78.2% |
| intracellular.ralstonia_like | 2,948 | 2,696 | 2,221 | 2,062 | 75.3% |
| intracellular.bacillus_like | 60 | 53 | 5 | 5 | 8.3% |
| extracellular.general | 596 | 500 | 404 | 365 | 67.8% |
| extracellular.lemoignei_like | 434 | 418 | 234 | 232 | 53.9% |
| **合计** | 8,731 | 7,051 | 6,532 | 5,531 | 74.8% |

> Bacillus 型验证率极低 (8.3%) 是因为大多数 c95 命中为弱同源、未通过 HMMer。

---

## 6. 古菌验证详情

| 指标 | 值 |
|---|---|
| 古菌基因组扫描 | 10,122 |
| DIAMOND 初始命中 | 2 (Halovenus sp. 34aa, UBA73 sp. 44aa) |
| 经过 ≥100aa 过滤后 | 0 |
| HMMer 验证后 | 0 |
| **结论** | **古菌中无高置信度 phaZ 基因** |

> 表述建议："在本次 GTDB R232 代表古菌基因组分析范围内，未发现高置信度 phaZ 基因。"

---

## 7. gRodon2 基因组集

| 分组 | 基因组数 | 来源 |
|---|---|---|
| phaZ+ (c95 级别匹配) | 4,394 | 含多拷贝基因组 |
| phaZ- (同属对照) | 4,394 | 同属但无 phaZ 命中 |
| 原始 manifest | 8,788 | 903 属 |
| 平衡后保留 | 8,692 | 剔除 53 失败 + 3 异常 |
| 失败原因 | 核糖体蛋白 HMM 命中不足 (<10) | — |

> gRodon2 使用 c95 级别基因组池（宽于验证集 5,531）以扩大配对范围。
