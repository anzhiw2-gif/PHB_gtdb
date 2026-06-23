# 可溯源可重复分析流程

本文档给出从环境、输入数据、脚本运行到结果检查的复现路线。服务器路径以 T141 当前项目为准，其他服务器可通过环境变量覆盖。

## 1. 运行环境

推荐服务器环境：

```text
Server: T141
Project path: /home/data/haoyu/PHB_gtdb
GTDB root: /home/data/haoyu/GTDB
GTDB release: R232
Python env: /home/data/haoyu/miniconda3/envs/phb_gtdb
R/gRodon env: /home/data/haoyu/miniconda3/envs/grodon2
```

创建 Python 环境：

```bash
conda env create -f environment.yml
conda activate phb_gtdb
```

设置路径：

```bash
export PHB_PROJECT_ROOT=/home/data/haoyu/PHB_gtdb
export GTDB_ROOT=/home/data/haoyu/GTDB
```

## 2. 输入数据

项目依赖 GTDB R232 代表基因组、taxonomy、metadata 和参考树：

```text
/home/data/haoyu/GTDB/gtdb_genomes_reps_r232/database
/home/data/haoyu/GTDB/taxonomy/bac120_taxonomy_r232.tsv
/home/data/haoyu/GTDB/taxonomy/ar53_taxonomy_r232.tsv
/home/data/haoyu/GTDB/metadata/bac120_metadata_r232.tsv.gz
/home/data/haoyu/GTDB/GTDB_tree/bac120_r232.tree
```

这些文件体积大，不进入 GitHub。仓库中保留的是脚本、轻量作图数据、文档和最终论文图。

## 3. 核心 PhaZ 分析

一键运行主流程：

```bash
cd /home/data/haoyu/PHB_gtdb
bash main_pipeline.sh --threads 30
```

单步运行：

```bash
bash main_pipeline.sh --step 1 --threads 30
bash main_pipeline.sh --step 2 --threads 30
bash main_pipeline.sh --step 3 --threads 30
bash main_pipeline.sh --step 4 --threads 30
bash main_pipeline.sh --step 5 --threads 30
bash main_pipeline.sh --step 6 --threads 30
```

每一步原理和脚本对应关系见 [SCRIPT_INDEX.md](SCRIPT_INDEX.md)。

## 4. 结果检查

运行：

```bash
python scripts/check_results.py
```

关键预期结果：

| 文件 | 预期数量 |
|---|---:|
| `data/processed/phaz_proteins_all.fasta` | 8,768 |
| `data/processed/phaz_proteins_c95.fasta` | 8,731 |
| `data/processed/phaz_proteins_filtered.fasta` | 7,478 |
| `data/processed/phaz_proteins_validated.fasta` | 6,532 |
| `data/processed/phb_search_results.tsv` | 7,068 rows, `phaZ_count` sum = 16,486 |
| `data/processed/archaea_phb_search_results.tsv` | 2 rows, all removed by downstream filters |

完整 manifest 见 [../RUN_MANIFEST.md](../RUN_MANIFEST.md)。

## 5. Nature 风格图复现

仓库已上传五张图及轻量源数据：

```text
figures/nature/
figure_data/
```

复现命令：

```bash
python scripts/07_nature_figures.py
```

输出：

```text
figures/nature/figure1_workflow_funnel.{pdf,png,svg}
figures/nature/figure2_phylum_heatmap.{pdf,png,svg}
figures/nature/figure3_subtype_lipase.{pdf,png,svg}
figures/nature/figure4_genera_phylogeny.{pdf,png,svg}
figures/nature/figure5_grodon_growth_comparison_hmm_allmatched.{pdf,png,svg}
```

图题和图注见 [FIGURE_CAPTIONS.md](FIGURE_CAPTIONS.md)。

## 6. gRodon2 生长速率分析

安装并验证 gRodon2：

```bash
Rscript scripts/install_grodon2.R /home/data/haoyu/software/gRodon2
```

准备 Pfam 原核核糖体蛋白 HMM 小库：

```bash
python scripts/08_prepare_ribosomal_hmms.py
```

小规模验证：

```bash
python scripts/08_grodon_growth.py \
  --heg-method hmm \
  --threads 2 \
  --marker-threads 1 \
  --max-per-genus 3 \
  --pilot 12
```

正式全匹配运行：

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

监控命令：

```bash
pgrep -af '08_grodon_growth.py'
tail -40 results/logs/grodon_growth_hmm_allmatched.out
wc -l results/tables/grodon_growth_predictions_hmm_allmatched.tsv
```

完成后生成严格同属平衡表、失败基因组汇总、统计检验和 Figure 5：

```bash
python scripts/10_balance_grodon_by_genus.py --label hmm_allmatched
python scripts/11_grodon_growth_stats.py --label hmm_allmatched
```

设计说明：

- 只在含有 `phaZ+` 的属中寻找同属 `phaZ-` 对照。
- 对每个属进行数量平衡，避免一边样本量过大。
- `--max-per-genus 0` 表示不设每属上限，使用所有可平衡配对的基因组。
- 全匹配 manifest 为 8,788 个基因组，即 4,394 个 `phaZ+` 和 4,394 个同属 `phaZ-`，覆盖 903 个属。
- gRodon2 完成后，8,735 个基因组 `status=ok`，53 个失败；其中 3 个 `status=ok` 记录缺失有效生长速率。
- 严格平衡后保留 8,692 个基因组，即 4,346 个 `phaZ+` 和 4,346 个 `phaZ-`，覆盖 899 个属。
- 属水平检验未发现显著全局差异；详细结果见 [GRODON2_FINAL_STATS.md](GRODON2_FINAL_STATS.md)。

## 7. 可追溯记录

建议每次正式运行后记录：

```text
日期
服务器
Git commit
conda env
命令
日志文件
核心输出表
结果校验状态
```

本次 Codex 协作整理见 [CODEX_WORKLOG.md](CODEX_WORKLOG.md)。
