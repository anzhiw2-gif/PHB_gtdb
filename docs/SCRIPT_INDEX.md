# 脚本索引与每一步原理

本文档用于把 `scripts/` 中的脚本和项目每一步分析逻辑对应起来，方便复现、检查和汇报。

## 总体原则

本项目的核心思想是先用高召回策略找到候选 PhaZ，再通过多重证据提高可信度，最后按亚型独立分析。原因是 PhaZ 属于 alpha/beta hydrolase 相关的功能多样酶家族，不同亚型之间序列差异较大，不适合不分型直接混合建树。

```text
候选搜索 -> 序列提取 -> 去冗余和过滤 -> HMM/DIAMOND 验证
        -> 亚型划分 -> 独立比对建树 -> 统计可视化 -> 生长速率扩展分析
```

## 核心脚本表

| 脚本 | 所属步骤 | 作用 | 主要原理 | 主要输出 |
|---|---|---|---|---|
| `config.py` | 全局配置 | 集中管理项目路径、GTDB 路径和参数 | 通过环境变量覆盖服务器路径，保证脚本可迁移 | 无直接输出 |
| `utils.py` | 工具函数 | FASTA 读写、命令运行、日志辅助等 | 减少重复代码，统一文件处理方式 | 无直接输出 |
| `01_phb_search.py` | Step 1 | 搜索细菌基因组中的 PhaZ 候选 | Pyrodigal 预测蛋白，DIAMOND 对 14 条验证 PhaZ 参考序列搜索同源 | `data/processed/phb_search_results.tsv` |
| `01b_archaea_search.py` | Step 1b | 搜索古菌基因组 | 同样使用 Pyrodigal + DIAMOND，用于检验古菌是否存在可信 PhaZ | `data/processed/archaea_phb_search_results.tsv` |
| `02_extract_sequences.py` | Step 2 | 提取 PhaZ 候选蛋白并过滤 | 根据 Step 1 命中坐标/ID 提取蛋白序列，进行 CD-HIT c95 去冗余和长度过滤 | `phaz_proteins_all.fasta`, `phaz_proteins_c95.fasta`, `phaz_proteins_filtered.fasta` |
| `03_msa.py` | Step 3 | 多序列比对 | 按亚型分别用 MAFFT 比对，再用 trimAl 修剪低质量区域 | `*_aligned.fasta`, `*_trim.fasta` |
| `04_phylogeny.py` | Step 4 | 系统发育树构建 | 小亚型用 IQ-TREE + ModelFinder + bootstrap，大亚型用 FastTree 解决计算量问题 | `*.treefile` |
| `05_annotation.py` | Step 5 | 功能与分类统计 | 整合 GTDB taxonomy、PhaZ 亚型、lipase box 和分布统计 | `results/tables/` |
| `06_visualization.py` | Step 6 | 常规可视化 | 生成长度分布、门水平分布、亚型组成、Top 属等图 | `results/figures/` |
| `07_nature_figures.py` | Figure | 生成四张 Nature 风格汇报图 | 读取 `figure_data/` 中的轻量源数据，输出 SVG/PDF/PNG | `figures/nature/figure*.{pdf,png,svg}` |
| `08_make_grodon_seed.R` | gRodon pilot | 从 gRodon2 示例基因组提取核糖体蛋白 seed | 用示例 CDS 名称中的 ribosomal protein 标记生成 seed，用于 DIAMOND pilot | `data/external/grodon/grodon_ribosomal_seed.*` |
| `08_prepare_ribosomal_hmms.py` | gRodon formal | 准备 Pfam 核糖体蛋白 HMM 小库 | 下载 Pfam metadata 和 HMM，筛选原核核糖体蛋白模型，抽取小型 HMM 库 | `data/external/grodon/ribosomal_pfam.hmm` |
| `08_run_grodon_one.R` | gRodon formal | 对单个基因组运行 gRodon2 | 读取 CDS 和 HEG ID，调用 `gRodon::predictGrowth` | 标准输出中的 `RESULT` 行 |
| `08_grodon_growth.py` | gRodon formal | 同属内 `phaZ+` vs `phaZ-` 最大生长速率预测 | Pyrodigal 预测 CDS，HMMER 标记核糖体 HEG，gRodon2 预测 doubling time，再按属汇总 | `results/tables/grodon_growth_*` |
| `09_monitor_grodon_progress.py` | gRodon monitor | 实时监控 gRodon2 后台任务 | 读取 manifest、prediction table 和 nohup log，计算进度、速度、ETA，并可输出 PNG 进度图 | 终端 dashboard, optional history TSV/PNG |
| `10_balance_grodon_by_genus.py` | gRodon statistics | 生成严格同属平衡比较表 | 只保留 `status=ok` 且有有效生长速率的记录；每个属取 `min(phaZ+, phaZ-)` 形成等量比较，并汇总失败基因组 | `grodon_growth_balanced_by_genus_*`, `grodon_failed_*` |
| `11_grodon_growth_stats.py` | gRodon statistics + figure | 统计检验并生成 Figure 5 | 以属为主分析单位，计算 `mean(phaZ+) - mean(phaZ-)`，进行 sign test、Wilcoxon、属内置换检验和 bootstrap CI | `grodon_growth_statistical_tests_*`, `figure5_grodon_growth_comparison_*` |
| `check_results.py` | QC | 检查关键结果是否齐全 | 按预期序列数、表格行数、tree 和 figure 文件存在性进行校验 | 终端 QC 报告 |

## Step 1: PhaZ 候选搜索

输入是 GTDB R232 代表基因组和 14 条经验证的 PhaZ 参考序列。脚本先用 Pyrodigal 对每个基因组预测蛋白编码基因，再用 DIAMOND 进行高通量同源搜索。这样做的好处是速度快，适合十万级基因组规模。

注意：Step 1 是候选搜索，不等于最终 PhaZ 集合。DIAMOND 初筛会保留一定边界命中，后续需要长度过滤、HMM 验证和高置信度标准进一步筛选。

## Step 2: 提取、去冗余和过滤

提取所有候选 PhaZ 蛋白后，使用 CD-HIT c95 去除近完全重复序列。c95 的目的不是压缩功能多样性，而是去掉几乎完全相同的重复序列。项目中不使用 c70/c60/c50 作为系统发育证据，因为这些阈值会把胞外、胞内等不同亚家族强制合并。

长度过滤用于去除明显不完整的短片段。古菌中的 2 条初始命中也在这里被识别为短片段，未进入最终验证集。

## Step 3-4: 亚型独立比对与建树

PhaZ 的亚型之间差异大，直接混合建树会让低同源区域主导树形。因此项目采用：

```text
先分型 -> 每个亚型独立 MAFFT 比对 -> trimAl 修剪 -> 单独建树
```

小规模亚型使用 IQ-TREE 和 ModelFinder，以获得更严格模型选择和 bootstrap 支持。`intracellular` 与 `ralstonia` 序列量大，使用 FastTree 作为计算可行的替代方案。

## Step 5-7: 统计、可视化和论文图

常规图保存在 `results/figures/`，Nature 风格汇报图保存在 `figures/nature/`。后者使用 `figure_data/` 中的轻量源数据，可在不携带完整 GTDB 和大中间文件的情况下复现图形。

## Step 8: gRodon2 生长速率扩展分析

专家建议比较同属内是否含有 `phaZ` 的菌株生长速率差异。为控制分类学和生态背景，脚本只在含 `phaZ+` 的属内寻找同属 `phaZ-` 对照，并进行数量平衡。

gRodon2 需要 CDS 和 highly expressed genes。项目使用 Pyrodigal 预测 CDS，用 Pfam 原核核糖体蛋白 HMM 识别核糖体蛋白，并把这些核糖体蛋白作为 HEG。该路线比仅使用 gRodon2 示例基因组的 DIAMOND seed 更适合作为正式分析。

输出包括：

- `grodon_growth_manifest_*.tsv`: 进入分析的同属配对清单
- `grodon_growth_predictions_*.tsv`: 每个基因组的 gRodon2 预测结果
- `grodon_growth_same_genus_summary_*.tsv`: 每个属的 `phaZ+` 与 `phaZ-` 均值差异
- `grodon_growth_balanced_by_genus_*.tsv`: 严格同属平衡后的基因组级比较表
- `grodon_growth_balanced_by_genus_summary_*.tsv`: 严格同属平衡后的属级均值差异表
- `grodon_growth_statistical_tests_*.tsv`: 属水平统计检验结果
- `figures/nature/figure5_grodon_growth_comparison_*.{pdf,png,svg}`: gRodon2 同属平衡分析图

## 结果校验

完成核心流程后运行：

```bash
python scripts/check_results.py
```

该脚本检查关键 FASTA 数量、搜索结果行数、树文件和图文件，便于快速发现漏跑或路径错误。
