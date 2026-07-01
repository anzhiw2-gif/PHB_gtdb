# PHB_gtdb 结果审计记录 2026-06-23

本审计基于两类来源：

1. GitHub 当前发布内容，即 `origin/main`。
2. T141 服务器结果目录 `/home/data/haoyu/PHB_gtdb`。

审计目标是检查项目关键结果是否完整、数量是否一致、文档是否存在过期描述，以及哪些环节需要在汇报中说明边界。

## 结论摘要

当前核心结果可以用于汇报和后续分析：

- PhaZ 主流程核心文件通过 `scripts/check_results.py` 检查。
- `phaz_proteins_all.fasta` 已恢复为 8,768 条序列，不再是 0 字节文件。
- 最终高置信度 PhaZ 集合为 6,532 条序列。
- gRodon2 全匹配任务已完成，平衡后得到 4,346 个 `phaZ+` 与 4,346 个 `phaZ-` 基因组，覆盖 899 个属。
- 属水平统计未发现 `phaZ+` 与 `phaZ-` 之间显著全局预测最大生长速率差异。

## 服务器检查结果

在 T141 运行：

```bash
cd /home/data/haoyu/PHB_gtdb
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/check_results.py
```

结果：

```text
All expected result files passed.
```

核心 FASTA 检查：

| 文件 | 序列数 | 状态 |
|---|---:|---|
| `phaz_proteins_all.fasta` | 8,768 | OK |
| `phaz_proteins_c95.fasta` | 8,731 | OK |
| `phaz_proteins_filtered.fasta` | 7,478 | OK |
| `phaz_proteins_validated.fasta` | 6,532 | OK |
| `phaz_bacillus_type.fasta` | 60 | OK |
| `phaz_extracellular.fasta` | 596 | OK |
| `phaz_extracellular_lemoignei.fasta` | 434 | OK |
| `phaz_intracellular.fasta` | 4,693 | OK |
| `phaz_ralstonia.fasta` | 2,948 | OK |

搜索结果检查：

| 文件 | 行数 | `phaZ_count` 总和 | 状态 |
|---|---:|---:|---|
| `phb_search_results.tsv` | 7,068 | 16,486 | OK |
| `archaea_phb_search_results.tsv` | 2 | 2 | OK |

## gRodon2 检查结果

核心结果文件存在：

```text
results/tables/grodon_growth_manifest_hmm_allmatched.tsv
results/tables/grodon_growth_predictions_hmm_allmatched.tsv
results/tables/grodon_growth_same_genus_summary_hmm_allmatched.tsv
results/tables/grodon_growth_balanced_by_genus_hmm_allmatched.tsv
results/tables/grodon_growth_balanced_by_genus_summary_hmm_allmatched.tsv
results/tables/grodon_growth_balanced_genus_selected_counts_hmm_allmatched.tsv
results/tables/grodon_growth_statistical_tests_hmm_allmatched.tsv
results/tables/grodon_growth_genus_effects_hmm_allmatched.tsv
results/tables/grodon_failed_genomes_hmm_allmatched.tsv
results/tables/grodon_failed_genera_summary_hmm_allmatched.tsv
```

最终统计：

| 指标 | 数量 |
|---|---:|
| gRodon2 输入基因组 | 8,788 |
| `status=ok` | 8,735 |
| failed | 53 |
| `status=ok` 但缺失有效生长速率 | 3 |
| 严格平衡后基因组 | 8,692 |
| 平衡后 `phaZ+` | 4,346 |
| 平衡后 `phaZ-` | 4,346 |
| 平衡后属数 | 899 |

统计检验：

| 检验 | 结果 |
|---|---:|
| 属水平平均差值 | -0.002734 h^-1 |
| 属水平中位差值 | +0.001813 h^-1 |
| Wilcoxon signed-rank P | 0.459 |
| Exact sign test P | 0.386 |
| 属内分层置换检验 P | 0.670 |

解释：严格同属平衡后，没有证据支持 `phaZ+` 与 `phaZ-` 基因组在预测最大生长速率上存在显著全局差异。

## 发现的问题和处理

### 1. README 中 gRodon2 状态已过期

旧描述写作“正在后台运行”，不符合当前结果。

处理：

- 已更新为“gRodon2 同属平衡比较已完成”。
- Figure 5 已并入原“已生成论文图”表格。
- 删除底部单独的英文 `Latest gRodon2 Growth-Rate Result` 补丁式小节。

### 2. 根目录存在早期临时监控脚本

根目录中原有：

```text
download_monitor.sh
extract_monitor.sh
monitor.sh
monitor_search.sh
```

这些脚本属于 GTDB 下载、解压和早期搜索监控阶段，不是当前主复现入口。

处理：

```text
scripts/legacy_monitoring/
```

并新增 `scripts/legacy_monitoring/README.md` 说明其历史用途。

### 3. Step 2 文件时间戳和下游文件时间戳不完全一致

服务器上 `phaz_proteins_all.fasta` 是 2026-06-23 重新提取得到，而 `phaz_proteins_c95.fasta`、`phaz_proteins_filtered.fasta`、`phaz_proteins_validated.fasta` 是早期正式运行结果。

判断：

- 重新提取后的 `phaz_proteins_all.fasta` 为 8,768 条，与历史正式结果一致。
- 下游 c95/filter/validated 文件数量通过 `check_results.py`，与预期完全一致。
- 因此当前汇报和统计结果可以使用。

建议：

- 汇报中无需展开“0 字节修复”过程。
- 如果未来要做方法论文级别完全严格复现，可从 Step 2 开始全链路重跑，使所有文件时间戳和输入链完全一致。

### 4. gRodon2 的失败记录不是整体错误

53 个 failed 主要来自核糖体 HMM 命中数不足，属于单个基因组不可预测，不代表全流程失败。

处理：

- 失败明细保存于 `grodon_failed_genomes_hmm_allmatched.tsv`。
- 正式统计只使用 `status=ok` 且有有效 `growth_rate_per_h` 的记录。

## 当前 GitHub 内容整理

本次整理后的 GitHub 结构更清晰：

```text
README.md                         项目入口和核心结果
RUN_MANIFEST.md                   运行位置、预期数量、最终结果账本
docs/PROJECT_TRACEABILITY.md      每一步脚本、代码原理、科学原理
docs/RESULTS_AUDIT_2026-06-23.md  本次结果审计记录
docs/GRODON2_FINAL_STATS.md       gRodon2 最终统计结果
docs/SCRIPT_INDEX.md              脚本索引
figures/nature/                   Figure 1-5
figures/nature/source_data/       Figure 1-5 source data
scripts/legacy_monitoring/        早期监控脚本归档
```

## 后续优先建议

1. 做 PhaZ 亚型分层的 gRodon2 分析，检查全局无差异是否由不同亚型方向抵消。
2. 对大样本属单独分析，如 `Bradyrhizobium`、`Rubrivivax`、`Polaromonas`、`Mesorhizobium`、`Polynucleobacter`。
3. 做系统发育相关性分析，避免把近缘基因组重复当作独立证据。
4. 如准备投稿级方法补充材料，建议全链路重跑 Step 2 之后的 c95、filter、validation 和 annotation，以获得完全一致的时间戳和 provenance。

