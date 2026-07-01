# Codex 协作记录摘要

说明：这里保存的是本次 VS Code/Codex 协作过程中形成的关键操作、决策和结果摘要，用于项目可追溯归档。它不是聊天窗口的逐字导出，但覆盖了当前项目整理、服务器运行、图形生成和 gRodon2 扩展分析的核心记录。

## 项目与服务器

```text
Local project: D:\codeX\PHB_gtdb
Archive copy requested: D:\PHB_gtdb
GitHub repo: https://github.com/anzhiw2-gif/PHB_gtdb
Server: T141
Server SSH: haoyu@10.16.1.141
Server project: /home/data/haoyu/PHB_gtdb
GTDB root: /home/data/haoyu/GTDB
```

## 已完成整理

1. 系统检查项目结构、脚本、结果和报告。
2. 整理汇报总结文档 `docs/PRESENTATION_SUMMARY.md`。
3. 安装并使用 `nature-figure` skill 生成四张 Nature 风格图。
4. 修改 Figure 1 左侧流程图，使流程方框完整且布局更美观。
5. 为 Figure 1-4 生成图题和图注。
6. 解释 lipase box 验证、胞外/胞内分型、`extracellular_lemoignei` 与一般胞外型差异。
7. 根据专家建议，制定三项扩展分析：
   - 同属内 `phaZ+` 与 `phaZ-` 最大生长速率比较
   - 亚型命名优化为“定位大类 + 谱系后缀”
   - 系统发育相关性和多样性分析

## 已生成图片

```text
figures/nature/figure1_workflow_funnel.pdf
figures/nature/figure1_workflow_funnel.png
figures/nature/figure1_workflow_funnel.svg
figures/nature/figure2_phylum_heatmap.pdf
figures/nature/figure2_phylum_heatmap.png
figures/nature/figure2_phylum_heatmap.svg
figures/nature/figure3_subtype_lipase.pdf
figures/nature/figure3_subtype_lipase.png
figures/nature/figure3_subtype_lipase.svg
figures/nature/figure4_genera_phylogeny.pdf
figures/nature/figure4_genera_phylogeny.png
figures/nature/figure4_genera_phylogeny.svg
```

## gRodon2 工作记录

已在 T141 下载并安装 gRodon2：

```text
gRodon2 source: /home/data/haoyu/software/gRodon2
R env: /home/data/haoyu/miniconda3/envs/grodon2
```

已验证 gRodon2 示例基因组可以运行 `predictGrowth()`。

先跑通 seed 版 pilot：

```text
Pyrodigal -> CDS/protein
DIAMOND -> gRodon2 示例核糖体 seed
gRodon2 predictGrowth
同属 phaZ+ vs phaZ-
```

seed 版 pilot 结果：

```text
12 rows, all status=ok
n_highly_expressed: 15-38
growth_rate_per_h: 0.0758-0.3449
```

随后升级为正式 HMM 路线：

```text
Pyrodigal -> CDS/protein
Pfam ribosomal protein HMM + HMMER -> HEG
gRodon2 predictGrowth
同属 phaZ+ vs phaZ-
```

严格 Pfam HMM pilot：

```text
12/12 genomes succeeded
n_highly_expressed: 35-45
mean growth_rate_per_h: 0.1787
```

全匹配 manifest：

```text
8,788 genomes
4,394 phaZ+
4,394 same-genus phaZ-
903 genera
```

正式后台任务已启动：

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

## 本次归档新增文件

```text
README.md
docs/SCRIPT_INDEX.md
docs/REPRODUCIBLE_WORKFLOW.md
docs/FIGURE_CAPTIONS.md
docs/GROWTH_RATE_ANALYSIS.md
docs/CODEX_WORKLOG.md
scripts/08_prepare_ribosomal_hmms.py
scripts/08_grodon_growth.py
scripts/08_run_grodon_one.R
scripts/08_make_grodon_seed.R
scripts/install_grodon2.R
scripts/test_grodon2.R
```

## 需要继续跟进

1. 等待 T141 上 `grodon_growth_hmm_allmatched` 完成。
2. 完成后生成最终统计表和 Figure 5 候选图。
3. 对 `phaZ+` 与 `phaZ-` 生长速率差异做属内统计检验。
4. 结合 GTDB 树或 PhaZ 树做系统发育相关性分析。
5. 若需要更严格 HEG 验证，可加入 GTDB-Tk/TIGRFAM marker 作为交叉验证。
