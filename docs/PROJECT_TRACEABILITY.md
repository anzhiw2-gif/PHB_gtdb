# PHB_gtdb 项目溯源总表

本文档用于回答三个问题：

1. 每一步由哪个脚本完成？
2. 脚本在代码层面做了什么？
3. 这一步背后的科学原理是什么？

服务器正式运行路径：

```text
/home/data/haoyu/PHB_gtdb
```

GTDB R232 数据路径：

```text
/home/data/haoyu/GTDB
```

## 总体逻辑

项目分为两条主线：

```text
PhaZ 鉴定主线:
GTDB genomes
  -> 基因预测
  -> PhaZ 同源搜索
  -> 候选蛋白提取
  -> 去冗余、长度过滤、HMM/DIAMOND 验证
  -> 亚型划分
  -> 独立比对、建树、统计和论文图

生长速率扩展主线:
同属 phaZ+ / phaZ- 基因组
  -> CDS 预测
  -> 核糖体蛋白 HMM 识别 HEG
  -> gRodon2 最大生长速率预测
  -> 同属严格平衡
  -> 属水平统计检验和 Figure 5
```

## Step 0. 全局配置与工具函数

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/config.py`, `scripts/utils.py` |
| 输入 | 环境变量、项目路径、GTDB 路径 |
| 输出 | 无直接结果文件 |
| 代码原理 | `config.py` 统一定义路径、线程数、阈值和 GTDB 文件位置；`utils.py` 统一 FASTA 读写、命令执行、日志和 taxonomy 解析 |
| 科学原理 | 保证同一套阈值和数据路径贯穿所有分析，减少手动运行时的路径漂移和参数漂移 |
| 质控点 | 修改服务器路径时优先使用环境变量覆盖，不直接改脚本内部常量 |

## Step 1. 细菌 PhaZ 候选搜索

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/01_phb_search.py` |
| 输入 | GTDB R232 细菌代表基因组；14 条经验证的 PhaZ 参考序列 |
| 输出 | `data/processed/phb_search_results.tsv` |
| 当前结果 | 7,068 个含 PhaZ 候选的细菌基因组；`phaZ_count` 总和 16,486 |
| 代码原理 | 对每个基因组用 Pyrodigal 预测蛋白，再用 DIAMOND blastp 对 PhaZ 参考库搜索同源命中，统计每个基因组的候选数量 |
| 科学原理 | PhaZ 属于同源可追踪的 PHB/PHA depolymerase 相关蛋白，先用高召回同源搜索扩大候选集合，再在下游严格验证 |
| 质控点 | Step 1 是候选搜索，不等于最终高置信度 PhaZ 集合 |

## Step 1b. 古菌 PhaZ 候选搜索

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/01b_archaea_search.py` |
| 输入 | GTDB R232 古菌基因组；同一 PhaZ 参考库 |
| 输出 | `data/processed/archaea_phb_search_results.tsv` |
| 当前结果 | 2 个初始短片段命中；最终未进入高置信度 PhaZ 集合 |
| 代码原理 | 与细菌搜索相同，独立扫描古菌数据以避免漏掉跨域同源 |
| 科学原理 | 用独立古菌扫描检验 PhaZ 是否可能存在于古菌；后续长度和保守特征过滤用于排除短片段假阳性 |
| 质控点 | 汇报中应表述为“古菌未发现高置信度 PhaZ”，而不是“完全无任何初始命中” |

## Step 2. PhaZ 候选蛋白提取、去冗余和过滤

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/02_extract_sequences.py` |
| 输入 | Step 1/1b 搜索结果；GTDB 基因组；PhaZ DIAMOND 库 |
| 输出 | `phaz_proteins_all.fasta`, `phaz_proteins_c95.fasta`, `phaz_proteins_filtered.fasta` |
| 当前结果 | all=8,768；c95=8,731；filtered=7,478 |
| 代码原理 | 根据命中基因组重新预测蛋白并 DIAMOND 搜索，提取命中的 PhaZ 候选；CD-HIT c95 去除近完全重复；长度过滤去掉明显短片段 |
| 科学原理 | c95 去冗余保留主要功能多样性，同时减少近重复序列对后续统计和树的影响；长度过滤降低片段化假阳性 |
| 质控点 | `phaz_proteins_all.fasta` 曾发生 0 字节损坏，已用 GitHub 版本 Step 2 重新提取恢复为 8,768 条。下游 c95/filter/validated 文件数量与历史正式结果一致，服务器 `check_results.py` 已通过。严格重新分析时可从 Step 2 开始全链路重跑以保持时间戳完全一致 |

## Step 3. 多序列比对

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/03_msa.py` |
| 输入 | 已分型或待分型的 PhaZ FASTA |
| 输出 | `*_aligned.fasta`, `*_trim.fasta` |
| 当前结果 | bacillus/extracellular/extracellular_lemoignei/intracellular/ralstonia 均有 trim FASTA |
| 代码原理 | 对每个 PhaZ 亚型分别调用 MAFFT 比对，再用 trimAl 删除低质量比对区域 |
| 科学原理 | PhaZ 亚型之间差异较大，混合比对容易让低同源区域主导信号；分型后独立比对更适合系统发育分析 |
| 质控点 | 不建议把所有 PhaZ 混合做一棵总树作为主要结论 |

## Step 4. 系统发育树构建

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/04_phylogeny.py` |
| 输入 | Step 3 修剪后的比对文件 |
| 输出 | `*.treefile` |
| 当前结果 | 小亚型使用 IQ-TREE；大亚型使用 FastTree 备选树 |
| 代码原理 | 对小型亚型运行 IQ-TREE/ModelFinder/bootstrap；对大型亚型使用 FastTree 保证计算可行性 |
| 科学原理 | 系统发育树用于检验 PhaZ 亚型是否形成相对稳定分支，并辅助解释不同类群中的扩张模式 |
| 质控点 | FastTree 树适合概览和大规模可视化，关键小分支可用 IQ-TREE 做更严格验证 |

## Step 5. 注释、亚型划分和统计

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/05_annotation.py` |
| 输入 | 验证后的 PhaZ FASTA；GTDB taxonomy；参考序列匹配和保守 motif 信息 |
| 输出 | `results/tables/` 中的分类、亚型、lipase box 和分布统计表 |
| 当前结果 | 最终验证集 6,532 条；覆盖 15 个细菌门、1,135 个属 |
| 代码原理 | 整合序列 ID、taxonomy、亚型、参考匹配、lipase box motif，并输出汇总表 |
| 科学原理 | PhaZ 是多样酶家族，不能只用“是否命中”解释功能；需要结合亚型、保守基序和分类背景 |
| 质控点 | 推荐使用“定位大类 + 谱系后缀”的亚型命名，如 `intracellular.cupriavidus_like` |

## Step 6. 常规可视化

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/06_visualization.py` |
| 输入 | Step 5 统计表 |
| 输出 | `results/figures/` |
| 代码原理 | 用 matplotlib/seaborn 生成长度分布、门水平分布、Top 属、亚型比例、lipase box 等常规图 |
| 科学原理 | 常规图用于探索数据结构和检查异常，Nature 图用于汇报和论文表达 |
| 质控点 | 常规图不等同于最终论文图；最终汇报图见 `figures/nature/` |

## Step 7. Nature 风格论文图

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/07_nature_figures.py` |
| 输入 | `figure_data/` 和 `figures/nature/source_data/` 中的轻量源数据 |
| 输出 | Figure 1-4 的 PDF/PNG/SVG |
| 代码原理 | 将工作流、漏斗、门水平热图、亚型组成、lipase box、Top 属和树概览组织为多面板图 |
| 科学原理 | 每张图服务一个结论：筛选流程可靠、分布不均匀、亚型有差异、系统发育存在稳定分支 |
| 质控点 | 图的源数据已随仓库保存，便于汇报复现；大体积原始数据不进入 GitHub |

## Step 8. gRodon2 核糖体 HMM 准备

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/08_prepare_ribosomal_hmms.py` |
| 输入 | Pfam-A metadata 和 HMM 库 |
| 输出 | `data/external/grodon/ribosomal_pfam.hmm` |
| 代码原理 | 从 Pfam 中筛选原核核糖体蛋白模型，排除 40S/60S、线粒体、质体和真核相关模型，抽取小型 HMM 库 |
| 科学原理 | gRodon2 需要 highly expressed genes；细菌核糖体蛋白通常稳定高表达，适合作为 HEG 参考 |
| 质控点 | Pfam 小库比扫描完整 Pfam 成本低，也比仅用示例 seed 更适合正式分析 |

## Step 9. gRodon2 最大生长速率预测

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/08_grodon_growth.py`, `scripts/08_run_grodon_one.R` |
| 输入 | 同属匹配 manifest；GTDB 基因组；ribosomal HMM |
| 输出 | `grodon_growth_manifest_hmm_allmatched.tsv`, `grodon_growth_predictions_hmm_allmatched.tsv`, `grodon_growth_same_genus_summary_hmm_allmatched.tsv` |
| 当前结果 | 8,788 个基因组进入分析；8,735 个 `status=ok`；53 个失败 |
| 代码原理 | Python 端负责 Pyrodigal CDS/蛋白预测、HMMER 标记核糖体 HEG、并行调度；R 端调用 gRodon2 `predictGrowth` 输出 doubling time 和 growth rate |
| 科学原理 | gRodon2 基于 codon usage bias 估计最大生长潜力，用于比较基因组层面的生长策略差异 |
| 质控点 | 结果是“预测最大生长速率”，不是实测培养速率；失败多来自核糖体 HMM 命中数不足 |

## Step 10. 同属严格平衡

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/10_balance_grodon_by_genus.py` |
| 输入 | `grodon_growth_predictions_hmm_allmatched.tsv` |
| 输出 | `grodon_growth_balanced_by_genus_hmm_allmatched.tsv`, `grodon_growth_balanced_by_genus_summary_hmm_allmatched.tsv`, `grodon_failed_*` |
| 当前结果 | 4,346 个 `phaZ+` 与 4,346 个 `phaZ-`；899 个属 |
| 代码原理 | 只保留 `status=ok` 且有有效 `growth_rate_per_h` 的记录；每个属取 `min(n_phaZ+, n_phaZ-)` 形成严格等量比较 |
| 科学原理 | 控制属水平分类学背景，避免一侧样本数过大造成比较偏倚 |
| 质控点 | 4 个属因只剩单侧成功样本被排除；40 个成功记录因严格平衡被剔除 |

## Step 11. 统计检验和 Figure 5

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/11_grodon_growth_stats.py` |
| 输入 | Step 10 平衡表和属级汇总表 |
| 输出 | `grodon_growth_statistical_tests_hmm_allmatched.tsv`, `grodon_growth_genus_effects_hmm_allmatched.tsv`, Figure 5 |
| 当前结果 | Wilcoxon P=0.459；sign test P=0.386；属内置换 P=0.670 |
| 代码原理 | 以属为主分析单位，计算 `delta = mean(phaZ+) - mean(phaZ-)`，再做 sign test、Wilcoxon、属内置换和 bootstrap CI |
| 科学原理 | 属内平衡后再进行属水平检验，避免把同属多个近缘基因组当作完全独立样本 |
| 质控点 | 当前结论是没有显著全局方向性差异；后续更适合按亚型、具体属或生态背景分层 |

## Step 12. 结果检查

| 项目 | 内容 |
|---|---|
| 脚本 | `scripts/check_results.py` |
| 输入 | 服务器上的核心 FASTA、TSV、tree 和 figure 文件 |
| 输出 | 终端 QC 报告 |
| 当前结果 | 2026-06-23 在 T141 运行通过，所有预期核心结果文件 passed |
| 代码原理 | 按预期数量检查 FASTA 序列数、搜索 TSV 行数和 `phaZ_count` 总和、tree 文件和常规 figure 文件 |
| 科学原理 | 通过数量和文件完整性快速发现漏跑、路径错或 0 字节结果 |
| 质控点 | check_results 是结构性质控，不替代生物学解释或统计检验 |

