# PHB_gtdb — 管线评审与改进计划

**评审日期**: 2026-06-09

---

## 1. 基因搜索 (Step 1)

### ✅ 当前做法
- 14 条 NCBI 验证 PhaZ 参考序列 (3门8属)
- Pyrodigal v3.7 进程内基因预测
- DIAMOND v2.1 blastp: `-e 1e-10 --id 30 --query-cover 50`
- 30 线程并行, chunksize=1, 增量保存 + 断点续传

### 📋 改进计划

| # | 建议 | 优先级 | 状态 |
|---|------|:---:|:---:|
| 1 | **HMMer + Pfam 二次验证**: 用 PF10503 (PHB depolymerase) + TIGR01849 (胞内型) 对 DIAMOND hits 做 hmmscan 确认, 降低假阳性 | 🔴 高 | 📋 TODO |
| 2 | **关键位点检测**: 对 DIAMOND hits 检测 lipase box (G-X-S-X-G) 或 catalytic triad 保守性, 标注完整性 | 🔴 高 | 📋 TODO |
| 3 | **比对覆盖度过滤**: 增加 `--query-cover 60` (当前50), 对 hits 再加 `alignment_length/query_length` 阈值 | 🟡 中 | 📋 TODO |
| 4 | **负对照**: 随机选取 100 个已知不含 PhaZ 的基因组 + 随机打乱序列作为负对照, 评估误报率 | 🟡 中 | 📋 TODO |

---

## 2. CD-HIT 去冗余

### ✅ 当前做法
- c95: 仅去真冗余 (≥95% identity), 8,768→8,731 (99.6%)
- c70 用于概览热图 (仅可视化, 不做系统发育证据)
- 同亚型内部 c85 去噪 (ralstonia 2,948→2,860; intracellular 4,693→4,522)

### 📋 改进计划

| # | 建议 | 优先级 | 状态 |
|---|------|:---:|:---:|
| 5 | **多拷贝处理策略**: 对同一基因组有多个 PhaZ 拷贝的情况, 统计拷贝数分布, 保留所有拷贝 (当前做法) 并标注在树中 | 🟡 中 | 📋 TODO |
| 6 | **拷贝数分布统计**: 作为功能注释的一部分, 输出每个基因组各亚型 PhaZ 的拷贝数 | 🟡 中 | 📋 TODO |
| 7 | **c85 阈值的生物学依据**: 同亚型内部 c85 基于序列保守性 (胞内型平均 identity ~45%), 高于此阈值的合并仍然合理 | 🟢 低 | ✅ 已记录 |

---

## 3. 结构域分型

### ✅ 当前做法
- 5 亚型: extracellular / extracellular_lemoignei / intracellular / ralstonia / bacillus_type
- 分类依据: DIAMOND 最佳匹配参考序列 + Pfam 结构域特征
- 各亚型独立建树 (MAFFT → trimAl → IQ-TREE)

### 📋 改进计划

| # | 建议 | 优先级 | 状态 |
|---|------|:---:|:---:|
| 8 | **SignalP/Phobius 验证胞外型**: 对所有 extracellular/extracellular_lemoignei 序列预测信号肽, 验证分泌特征 | 🔴 高 | 📋 TODO |
| 9 | **Lipase box Logo 图**: 对各亚型提取 lipase box (G-X-S-X-G) 区域做序列保守性分析, 生成 WebLogo 图, 用于结果解释 | 🟡 中 | 📋 TODO |
| 10 | **AlphaFold2 结构预测**: 对 bacillus_type (60条) 和古菌新发现的 Halovenus/UBA73 序列做结构预测, 验证是否可功能化为 PhaZ | 🟢 低 | 📋 TODO |
| 11 | **古菌序列 reciprocal BLASTP**: 对 2 个古菌命中做 NCBI nr 回检, 排除 α/β 水解酶交叉匹配 | 🔴 高 | 📋 TODO |

---

## 4. 多序列比对与建树

### ✅ 当前做法
- MAFFT `--auto` 自动选择最优方法
- trimAl `-gappyout` 自动去除 gap-rich 区域
- IQ-TREE `-m MFP` (ModelFinder Plus) 自动模型选择
- `-B 1000` 标准 bootstrap
- 小树: bacillus_type (8线程) / extracellular_lemoignei (10) / extracellular (10)
- 大树: ralstonia (30) / intracellular (30)

### 📋 改进计划

| # | 建议 | 优先级 | 状态 |
|---|------|:---:|:---:|
| 12 | **模型选择文档化**: README 中提及 IQ-TREE 使用 `-m MFP` (ModelFinder Plus), 自动选择最优氨基酸替代模型 | 🟢 低 | ✅ 已实现 |
| 13 | **SH-aLRT 支持度**: IQ-TREE 默认输出 SH-aLRT + UFBoot, 树文件中已包含分支支持度 | 🟢 低 | ✅ 已实现 |
| 14 | **联合树 vs 单型树对比**: 构建所有亚型合并树, 检验分型合理性 (同亚型应形成单系群) | 🟡 中 | 📋 TODO |
| 15 | **基因树-物种树 reconciliation**: 对 intracellular 大组做 reconciliation 分析, 评估 HGT 事件 | 🟢 低 | 📋 TODO |

---

## 5. 功能注释

### 📋 改进计划

| # | 建议 | 优先级 | 状态 |
|---|------|:---:|:---:|
| 16 | **Pfam/InterProScan 功能域注释**: 对所有 8,731 条 c95 序列跑 InterProScan, 标注 Pfam 结构域组成 | 🔴 高 | 📋 TODO |
| 17 | **胞内外定位统计**: 结合 SignalP 结果统计胞外/胞内比例, 按门分组可视化 | 🟡 中 | 📋 TODO |
| 18 | **底物特异性预测**: 文献中有部分 PhaZ 底物特异性预测方法, 可作为探索性分析 | 🟢 低 | 📋 TODO |

---

## 6. 可视化

### 📋 改进计划

| # | 建议 | 优先级 | 状态 |
|---|------|:---:|:---:|
| 19 | **功能分布热图**: 不同 PhaZ 亚型在门/属层级的分布热图 (phyloseq 风格) | 🔴 高 | 📋 TODO |
| 20 | **重要菌株标注**: 对 Cupriavidus, Paucimonas (参考来源属), 古菌命中等在树中高亮标注 | 🟡 中 | 📋 TODO |
| 21 | **c70 概览热图 + 亚型标注**: c70 序列聚类的热图, 每个聚类标注对应亚型颜色 | 🟡 中 | 📋 TODO |

---

## 优先级汇总

### 🔴 高优先级 (下一阶段执行)
1. HMMer Pfam 二次验证 (降低假阳性)
2. Lipase box / catalytic triad 关键位点检测
3. SignalP 验证胞外型信号肽
4. 古菌序列 reciprocal BLASTP 验证
5. InterProScan 功能域注释
6. 功能分布热图

### 🟡 中优先级
7. 负对照评估误报率
8. 多拷贝统计
9. Lipase box Logo 图
10. 联合树 vs 单型树对比
11. 胞内外定位统计

### 🟢 低优先级 (探索性)
12. AlphaFold2 结构预测
13. 基因树-物种树 reconciliation
14. 底物特异性预测
15. c85 阈值生物学依据 (已记录)

---

## 当前状态

| 步骤 | 状态 |
|---|---|
| Step 1: 搜索 | ✅ 完成 |
| Step 2: 提取 | ✅ 完成 |
| CD-HIT c95 | ✅ 完成 |
| 结构域分型 | ✅ 完成 |
| 3 小组 IQ-TREE | ✅ 完成 |
| 2 大组 IQ-TREE | 🔄 进行中 |
| HMMer 二次验证 | 📋 TODO |
| SignalP + Logo | 📋 TODO |
| InterProScan | 📋 TODO |
| 可视化 | 📋 TODO |

---

**最后更新**: 2026-06-09
