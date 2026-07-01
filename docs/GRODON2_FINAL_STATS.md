# gRodon2 同属平衡生长速率分析最终结果

本文件总结专家建议中“同属内比较 `phaZ+` 与 `phaZ-` 菌株最大生长速率”的正式分析结果。这里的生长速率来自 gRodon2 基于 codon usage bias 的最大生长速率预测，不等同于特定培养条件下的实测生长速率。

## 分析问题

核心问题是：在同一个属内，携带 `phaZ` 的基因组是否相对于不携带 `phaZ` 的同属基因组表现出不同的预测最大生长速率？

为了避免不同属之间生活史、GC 含量、生态位和系统发育背景差异造成混杂，本分析不直接跨所有基因组比较，而是在每个属内做 `phaZ+` 与 `phaZ-` 数量平衡后，再进行属水平统计检验。

## 数据账本

| 项目 | 数量 |
|---|---:|
| gRodon2 输入基因组 | 8,788 |
| gRodon2 失败基因组 | 53 |
| 失败涉及属 | 37 |
| `status=ok` 但缺失有效生长速率 | 3 |
| 可用于生长速率比较的成功记录 | 8,732 |
| 为严格同属平衡额外剔除的成功记录 | 40 |
| 最终平衡分析基因组 | 8,692 |
| 最终 `phaZ+` 基因组 | 4,346 |
| 最终 `phaZ-` 基因组 | 4,346 |
| 最终覆盖属 | 899 |

平衡规则为：

```text
每个属内选取数量 = min(n_phaZ_positive_ok, n_phaZ_negative_ok)
```

因此每个纳入分析的属都满足：

```text
n(phaZ+) = n(phaZ-)
```

## 失败原因

gRodon2 失败的主要原因是核糖体蛋白 HMM 命中数不足：

| 错误类型 | 数量 |
|---|---:|
| `too_few_ribosomal_hits:9` | 26 |
| `too_few_ribosomal_hits:8` | 12 |
| `too_few_ribosomal_hits:7` | 8 |
| `too_few_ribosomal_hits:6` | 5 |
| `too_few_ribosomal_hits:4` | 1 |
| `hmmsearch_failed` | 1 |

失败最多的属包括 `Polynucleobacter`、`Limnohabitans_A`、`Polaromonas`、`Bradyrhizobium` 和 `Allorhodoferax`。

## 统计设计

主分析单位为“属”，而不是单个基因组。对每个属计算：

```text
delta = mean(growth_rate_per_h of phaZ+) - mean(growth_rate_per_h of phaZ-)
```

然后对 899 个属的 `delta` 做总体检验。

采用的统计检验：

1. Wilcoxon signed-rank test：检验属水平差值是否系统性偏离 0。
2. Exact sign test：只看方向，检验 `delta > 0` 与 `delta < 0` 的属数是否偏离 1:1。
3. 属内分层置换检验：在每个属内部随机打乱 `phaZ+`/`phaZ-` 标签，评估观察到的平均差值是否超出随机期望。
4. Bootstrap 95% CI：对属水平均值、中位数和按样本数加权均值估计置信区间。

## 主要结果

| 指标 | 结果 |
|---|---:|
| 属水平平均差值 | -0.002734 h^-1 |
| 平均差值 bootstrap 95% CI | -0.015509 到 0.009320 h^-1 |
| 属水平中位差值 | +0.001813 h^-1 |
| 中位差值 bootstrap 95% CI | -0.002333 到 0.005652 h^-1 |
| 按每属配对数加权平均差值 | +0.008233 h^-1 |
| 加权平均 bootstrap 95% CI | -0.003591 到 0.019790 h^-1 |
| `phaZ+` 更高的属 | 463 |
| `phaZ-` 更高的属 | 436 |
| Exact sign test P | 0.386 |
| Wilcoxon signed-rank P | 0.459 |
| Wilcoxon effect r | 0.0247 |
| 属内分层置换检验 P，非加权平均 | 0.670 |
| 属内分层置换检验 P，加权平均 | 0.231 |

描述性基因组水平结果：

| 指标 | `phaZ+` | `phaZ-` |
|---|---:|---:|
| 平均 `growth_rate_per_h` | 0.391742 | 0.383509 |
| 中位 `growth_rate_per_h` | 0.239962 | 0.231760 |

基因组水平均值仅作为描述性结果，不作为主检验依据，因为同属内基因组不是完全独立样本。

## 结论

在严格同属平衡后，当前数据没有支持 `phaZ+` 基因组在总体上具有显著更高或更低预测最大生长速率的证据。属水平差值的中位数接近 0，Wilcoxon、sign test 和属内分层置换检验均不显著。

更合适的表述是：

> 在 GTDB R232 的同属平衡比较中，`phaZ` 携带状态与 gRodon2 预测最大生长速率之间未观察到一致的全局方向性差异；这提示 `phaZ` 更可能与特定类群、生态背景或 PhaZ 亚型相关，而不是简单对应整体更快或更慢的最大生长潜力。

## 生成文件

统计脚本：

```text
scripts/10_balance_grodon_by_genus.py
scripts/11_grodon_growth_stats.py
```

核心表格：

```text
results/tables/grodon_growth_balanced_by_genus_hmm_allmatched.tsv
results/tables/grodon_growth_balanced_by_genus_summary_hmm_allmatched.tsv
results/tables/grodon_growth_balanced_genus_selected_counts_hmm_allmatched.tsv
results/tables/grodon_growth_genus_effects_hmm_allmatched.tsv
results/tables/grodon_growth_statistical_tests_hmm_allmatched.tsv
results/tables/grodon_failed_genomes_hmm_allmatched.tsv
results/tables/grodon_failed_genera_summary_hmm_allmatched.tsv
```

论文图：

```text
figures/nature/figure5_grodon_growth_comparison_hmm_allmatched.svg
figures/nature/figure5_grodon_growth_comparison_hmm_allmatched.pdf
figures/nature/figure5_grodon_growth_comparison_hmm_allmatched.png
```

Figure 5 source data：

```text
figures/nature/source_data/figure5_grodon_balanced_genomes_hmm_allmatched.tsv
figures/nature/source_data/figure5_grodon_genus_effects_hmm_allmatched.tsv
figures/nature/source_data/figure5_grodon_statistical_tests_hmm_allmatched.tsv
```

## 后续建议

1. 按 PhaZ 亚型分层分析：重点比较 `intracellular.*`、`extracellular.*` 是否有不同生长速率模式。
2. 对大样本属单独建模：如 `Bradyrhizobium`、`Rubrivivax`、`Polaromonas`、`Mesorhizobium`、`Polynucleobacter`。
3. 结合系统发育相关性：使用系统发育距离或属/科水平随机效应，避免近缘基因组重复放大信号。
4. 将 gRodon2 结果与基因组 GC、CDS 数、生态来源或 PHB/PHA 相关代谢基因簇共同分析。
5. 汇报时避免说“实测生长速率”，建议使用“预测最大生长速率”或“最大生长潜力预测”。
