# PHB_gtdb — 硕士论文 Chapter 2 补充分析设计文档

**日期**: 2026-07-02
**上下文**: 项目所有核心分析已完成，此文档规划硕士论文 Chapter 2 的补充分析、论文结构、时间线
**运行平台**: T141 (10.16.1.141), `/home/data/haoyu/PHB_gtdb`

---

## 1. 硕士论文整体结构

```
论文标题（建议）: PHB 塑料降解菌的基因组全景扫描、功能验证与环境分布

Chapter 1  引言与文献综述
Chapter 2  基于 GTDB 全景扫描的 PhaZ 基因鉴定与进化分析  ← 本设计文档聚焦
Chapter 3  PHB 降解菌的分离培养与酶活性验证             ← 湿实验（不在本范围）
Chapter 4  环境宏基因组中 PhaZ 的多样性与分布特征        ← 宏基因组（不在本范围）
Chapter 5  综合讨论与结论
```

**三章逻辑链**: Chapter 2（全景→告诉你在哪里找）→ Chapter 3（实验室验证）→ Chapter 4（自然界验证）

---

## 2. Chapter 2 内容分解

| 小节 | 标题 | 状态 | 预计字数 |
|---|---|---|---|
| 2.1 | PhaZ 候选序列的系统搜索与验证 | ✅ 完成 | ~1500 |
| 2.2 | PhaZ 亚型分类与功能特征 | ✅ 完成 | ~2000 |
| 2.3 | PhaZ 系统发育与进化分析 | 🔄 70% — 需补充 T1+T2 | ~2500 |
| 2.4 | PhaZ 蛋白结构与功能预测 | ❌ 全新增 — T3+T4+T5 | ~2000 |
| 2.5 | gRodon2 预测生长速率分析 | ✅ 完成 | ~1500 |
| 2.6 | 宏基因组搜索种子序列准备 | ❌ 全新增 — T6 | ~800 |

---

## 3. 需要补充的 7 个任务

### T1 — phaZ 拷贝数分布统计

- **输入**: `data/processed/phb_search_results.tsv`
- **方法**: 统计每基因组 phaZ_count 分布；按门/属汇总平均拷贝数；高拷贝属文献比对
- **Skill**: 无需特殊 skill（纯 Python+pandas）
- **预计耗时**: 0.5h
- **输出**: 拷贝数分布直方图 + 门/属级汇总表 → 用于 Fig 2.6A

### T2 — PhaC 搜索与 phaZ-phaC 共进化分析

- **输入**: ~10 条 PhaC 参考序列（从文献收集 Class I-IV）；GTDB 基因组
- **步骤**:
  1. `protein-sequence-similarity-search`（MMseqs2）搜索 GTDB 中 PhaC 同源
  2. 筛选 PhaC 蛋白序列（E≤1e-10, ≥30% identity）
  3. 交叉 phaZ 和 phaC 的基因组 ID，做共分布矩阵
  4. `phylogenetics` skill 对 PhaC 序列建树
  5. 比较 PhaZ 树（已有）与 PhaC 树的属级别拓扑一致性
- **Skill**: `protein-sequence-similarity-search`, `phylogenetics`
- **预计耗时**: 3-4h
- **输出**: PhaC 树 + phaZ-phaC 共分布热图 → 用于 Fig 2.6B

### T3 — 5 亚型代表序列 3D 结构预测 【暂缓】

- **状态**: ⏸️ 暂缓。T141 无 GPU，ESMFold CPU 版可在需要时补充。
- **决定**: 优先完成其他 6 个任务，硕士论文核心故事不需要 3D 结构。

### T4 — 催化域保守性深度分析

- **输入**: 各亚型已有 MAFFT 比对 + trimAl 修剪文件（`data/processed/*_trim.fasta`）
- **方法**: `protein-sequence-msa`（Clustal Omega 重比对催化域区域 ~100aa）；生成序列 logo
- **Skill**: `protein-sequence-msa`, `nature-figure`
- **预计耗时**: 1h
- **输出**: 5 亚型催化域序列 logo 对比图 → 用于 Fig 2.7B

### T5 — 信号肽与分泌机制预测

- **输入**: `phaz_proteins_validated.fasta`（6,532 条）
- **方法**: SignalP 6.0（`conda install -c bioconda signalp6`）
- **分析**: 胞外型 vs 胞内型信号肽检出率；lemoignei 型分泌信号强度
- **Skill**: 无直接 skill（命令行工具）
- **预计耗时**: 0.5h
- **输出**: 信号肽统计表 + 分型对比图 → 用于 Fig 2.8A

### T6 — 宏基因组种子序列选拔

- **输入**: 6,532 条 PhaZ → CD-HIT c70 压缩 → 每亚型选代表性序列
- **方法**: `protein-sequence-similarity-search` 或手动 CD-HIT；每亚型选 5-10 条高多样性序列
- **Skill**: `protein-sequence-similarity-search`
- **预计耗时**: 1h
- **输出**: 精简 PhaZ 参考 FASTA（~50-100 条）→ 传给 Chapter 4 → 用于 Fig 2.8B

### T7 — 汇总重绘 Nature 图

- **方法**: `nature-figure` skill，将所有新分析结果更新到 Figure 中
- **Skill**: `nature-figure`
- **预计耗时**: 1h
- **输出**: Fig 2.6, 2.7, 2.8（各 PDF+PNG+SVG）

---

## 4. Figure 规划

| 图号 | 内容 | 来源 | 状态 |
|---|---|---|---|
| Fig 2.1 | 工作流 + 筛选漏斗 | 已有 Figure 1 | ✅ |
| Fig 2.2 | 门水平分布 + 亚型热图 | 已有 Figure 2 | ✅ |
| Fig 2.3 | 亚型组成 + Lipase box | 已有 Figure 3 | ✅ |
| Fig 2.4 | 属分布 + 系统发育树概览 | 已有 Figure 4 | ✅ |
| Fig 2.5 | gRodon2 生长速率比较 | 已有 Figure 5 | ✅ |
| Fig 2.6 | **拷贝数分布 + PhaZ-PhaC 共进化** | T1 + T2 | 🔄 新 |
| Fig 2.7 | **5亚型 3D 结构对比 + 催化域 Logo** | T3 + T4 | ❌ 新 |
| Fig 2.8 | **信号肽预测 + 宏基因组种子序列** | T5 + T6 | ❌ 新 |

---

## 5. T141 执行路线

### 第一轮（Day 1）
```
T1 拷贝数分布 → 0.5h
T5 SignalP → 0.5h
T4 催化域 Logo（可与 T1/T5 并行）→ 1h
```
全部无需 GPU，即跑即出。

### 第二轮（Day 2-3）
```
T3 AlphaFold2 → 后台 GPU job（1-2 天）
T2 PhaC 搜索 + 共进化 → 前台跑（~4h）
T6 种子序列 → 前台跑（1h）
```
T3 与 T2/T6 并行。

### 第三轮（Day 4）
```
T7 汇总出图 → 1h
开始 Chapter 2 写作
```

---

## 6. 论文写作阶段

| 阶段 | 内容 | 时间 |
|---|---|---|
| Phase 1 | Chapter 2 初稿（基于已完成数据） | Day 1-3 |
| Phase 2 | 补充分析结果写入 2.3-2.4-2.6 | Day 4-5 |
| Phase 3 | Chapter 2 完整稿 + 所有 Figure | Day 6-7 |
| Phase 4 | Chapter 1/3/4/5（配合导师） | Aug-Dec 2026 |
| Phase 5 | 全文统稿、参考文献、格式 | Jan-Feb 2027 |
| Deadline | 提交毕业论文 | 2027年3月 |

---

## 7. 依赖与风险

| 风险 | 影响 | 缓解措施 |
|---|---|---|
| T141 无 GPU | T3 无法跑 AlphaFold2 | 用 ESMFold CPU 版或 ColabFold 在线版 |
| PhaC 参考序列不足 | T2 范围受限 | 从 UniProt 和文献扩展参考集 |
| SignalP 6.0 安装困难 | T5 无法进行 | 备用 DeepSig 或 Phobius |
| 时间不够 | T3 可降级为仅 1-2 个亚型结构 | 选择差异最大的胞内/胞外代表 |

---

## 8. 技能使用映射

| Skill | 任务 | 用途 |
|---|---|---|
| `phylogenetics` | T2 | PhaC 建树 |
| `protein-sequence-similarity-search` | T2, T6 | MMseqs2 搜索 PhaC / CD-HIT 压缩 |
| `protein-sequence-msa` | T4 | 催化域 Clustal Omega 比对 |
| `nature-figure` | T4, T7 | 序列 Logo 图 + 最终汇总 |
| `statistics-verifier` | T1 | 拷贝数统计验证 |
| `bioinformatics-fundamentals` | T1-T6 | 数据格式检查 |
