# gRodon2 生长速率扩展分析

本分析对应专家建议：比较同属内携带 `phaZ` 与不携带 `phaZ` 的菌株是否存在最大生长速率差异。

## 科学问题

PhaZ 参与 PHB/PHA 降解和碳储存利用。若携带 `phaZ` 的菌株在生态策略上与同属不携带菌株不同，可能表现为预测最大生长速率差异。由于不同属之间生活史和基因组特征差异很大，本项目采用同属内比较，而不是跨所有细菌直接比较。

## 分析设计

```text
只选择含有 phaZ+ 基因组的属
  -> 在同属内寻找 phaZ- 对照基因组
  -> 每个属内进行 phaZ+ 和 phaZ- 数量平衡
  -> 对所有入选基因组预测最大生长速率
  -> 按属汇总 phaZ+ 与 phaZ- 差异
```

当前全匹配 manifest：

```text
8,788 个基因组
4,394 个 phaZ+
4,394 个同属 phaZ-
903 个属
```

## 为什么使用核糖体蛋白作为 HEG

gRodon2 需要两类输入：

1. 每个基因组的 CDS nucleotide sequences
2. highly expressed genes, 即高表达基因集合

在细菌中，核糖体蛋白基因通常稳定高表达，且跨物种高度保守，是 gRodon/codon usage bias 分析中常用的 HEG 参照。本项目因此使用核糖体蛋白作为 HEG。

## 为什么使用 Pfam ribosomal protein HMM

当前正式路线使用 Pfam 中的原核核糖体蛋白模型：

```text
Pfam-A metadata
  -> 筛选 ribosomal protein 模型
  -> 排除 40S/60S、线粒体、质体、修饰酶和真核相关模型
  -> 抽取小型 ribosomal_pfam.hmm
  -> HMMER 搜索每个基因组的预测蛋白
```

这样做的优点：

- HMM 模型来源公开且可复现。
- HMMER 对远缘同源比单纯 DIAMOND seed 更稳健。
- 只扫描核糖体蛋白小库，避免扫描完整 Pfam 的高成本。
- 不依赖完整 GTDB-Tk reference data。

GTDB-Tk/TIGRFAM marker 可作为后续交叉验证来源，但不是 gRodon2 的必要前提。

## 技术路线

```text
GTDB genome
  -> Pyrodigal 预测 CDS 和蛋白
  -> HMMER + Pfam ribosomal protein HMM 识别核糖体蛋白
  -> 核糖体蛋白 ID 作为 HEG
  -> gRodon2 predictGrowth
  -> 输出 doubling time 和 growth rate
```

对应脚本：

- `scripts/08_prepare_ribosomal_hmms.py`
- `scripts/08_grodon_growth.py`
- `scripts/08_run_grodon_one.R`

## 已完成 pilot

严格 Pfam HMM 版 pilot：

```text
12/12 genomes succeeded
n_highly_expressed: 35-45
mean n_highly_expressed: 41.3
growth_rate_per_h: 0.0746-0.5207
mean growth_rate_per_h: 0.1787
```

pilot 同属结果：

| Genus | n phaZ+ | n phaZ- | mean growth phaZ+ | mean growth phaZ- | delta |
|---|---:|---:|---:|---:|---:|
| `0-14-3-00-62-12` | 3 | 3 | 0.1532 | 0.1122 | +0.0409 |
| `17J80-11` | 1 | 1 | 0.4131 | 0.5207 | -0.1075 |
| `2-02-FULL-66-14` | 2 | 2 | 0.1151 | 0.0923 | +0.0228 |

pilot 只用于验证流程可行性，不能作为最终生物学结论。

## 正式运行

T141 后台任务：

```bash
nohup python scripts/08_grodon_growth.py \
  --heg-method hmm \
  --threads 8 \
  --marker-threads 1 \
  --max-per-genus 0 \
  --resume \
  --output-label hmm_allmatched \
  > results/logs/grodon_growth_hmm_allmatched.out 2>&1 &
```

输出文件：

```text
results/tables/grodon_growth_manifest_hmm_allmatched.tsv
results/tables/grodon_growth_predictions_hmm_allmatched.tsv
results/tables/grodon_growth_same_genus_summary_hmm_allmatched.tsv
results/logs/grodon_growth_hmm_allmatched.out
```

正式运行已完成。最终状态为：

```text
8,788 个基因组进入 gRodon2 全匹配分析
8,735 个基因组 status=ok
53 个基因组 failed
3 个 status=ok 记录缺失有效 growth_rate_per_h
8,692 个基因组进入严格同属平衡统计
4,346 个 phaZ+
4,346 个 phaZ-
899 个属
```

平衡与统计检验命令：

```bash
python scripts/10_balance_grodon_by_genus.py --label hmm_allmatched
python scripts/11_grodon_growth_stats.py --label hmm_allmatched
```

最终统计结果见 [GRODON2_FINAL_STATS.md](GRODON2_FINAL_STATS.md)。核心结论是：严格同属平衡后，未检测到 `phaZ+` 与 `phaZ-` 基因组之间显著的全局预测最大生长速率差异。

监控：

```bash
pgrep -af '08_grodon_growth.py'
tail -40 results/logs/grodon_growth_hmm_allmatched.out
wc -l results/tables/grodon_growth_predictions_hmm_allmatched.tsv
```

实时可视化监控：

```bash
python scripts/09_monitor_grodon_progress.py \
  --interval 60 \
  --history results/tables/grodon_growth_monitor_history.tsv \
  --plot results/figures/grodon_growth_progress.png
```

该脚本会显示终端进度条、完成比例、`ok/failed/pending` 数量、最近日志、处理速度和预计完成时间。`--plot` 会持续刷新一张 PNG 进度图，`--history` 会保存每次刷新记录，便于之后追踪运行速度变化。

## 后续统计建议

正式结果完成后建议进行：

1. 属内配对/平衡比较，统计 `delta = mean(phaZ+) - mean(phaZ-)`。
2. 对每个属计算效应量和置信区间。
3. 使用置换检验或分层模型评估总体方向。
4. 分亚型比较 `intracellular.*` 和 `extracellular.*` 是否对应不同生长策略。
5. 与系统发育距离结合，检查 `phaZ` 分布是否具有系统发育相关性。

## 解释边界

gRodon2 输出是基于 codon usage bias 的预测最大生长速率，不等同于某一培养条件下实测生长速率。因此最终表述应为“预测最大生长潜力”或“最大生长速率预测”，不应表述为实际培养速率。
