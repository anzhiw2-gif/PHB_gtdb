# 补充分析采样策略记录

**日期**: 2026-07-02 | **目的**: 记录所有补充分析中的抽样决策与方法

---

## 1. T3: 3D 结构预测（ESM Atlas API）

| 亚型 | 验证总数 | 预测数 | 采样策略 | 截取长度 |
|---|---|---|---|---|
| intracellular.cupriavidus_like | 3,668 | 5 (1+4) | 选最长序列 | 300aa |
| intracellular.ralstonia_like | 2,221 | 5 (1+4) | 选最长序列 | 300aa |
| intracellular.bacillus_like | 5 | 5 (1+4) | 全部预测 | 109-207aa 全长 |
| extracellular.general | 404 | 5 (1+4) | 选最长序列 | 300aa |
| extracellular.lemoignei_like | 234 | 5 (1+4) | 选最长序列 | 200-300aa |
| **合计** | 6,532 | **28** (23+5) | | |

**策略理由**:
- 每亚型选 5 条最长序列（最长 = 最完整蛋白序列），其中第一批（Phase 1）2 条覆盖胞内/胞外极端对比
- 截取 200-400aa 催化域：ESM Atlas API 限制 ~400aa；催化域通常在 N 端前 400aa
- Bacillus 型仅 5 条，全部预测（2 条 504 超时失败）
- 方法引用: ESMFold v1 (Lin et al., 2023, *Science*)

---

## 2. #14: 联合系统发育树（单系性检验）

| 亚型 | 修剪后总数 | 采样数 | 比例 | 策略 |
|---|---|---|---|---|
| intracellular.cupriavidus_like | 4,424 | 250 | 5.7% | 比例抽稀 |
| intracellular.ralstonia_like | 2,776 | 250 | 9.0% | 比例抽稀 |
| intracellular.bacillus_like | 57 | **57** | **100%** | 全部保留 |
| extracellular.general | 509 | 150 | 29.5% | 比例抽稀 |
| extracellular.lemoignei_like | 412 | 150 | 36.4% | 比例抽稀 |
| **合计** | 8,178 | **857** | 10.5% | |

**策略理由**:
- Bacillus 型仅 57 条，全部保留——抽样会导致其失去代表性
- 大亚型（>400 条）抽稀至 150-250 条以控制 MAFFT 计算时间在 ~30 分钟
- 中等亚型（400-500 条）保留 ~30-37% 以维持多样性
- 随机种子 42 确保可重复
- 方法: MAFFT `--auto` → trimAl `-automated1` → FastTree `-lg -gamma`

---

## 3. #9: Lipase box 序列 Logo 分析

| 亚型 | 比对序列数 | 含 G-X-S-X-G 数 | 比例 | 说明 |
|---|---|---|---|---|
| intracellular.cupriavidus_like | 4,424 | 24 | 0.5% | trimAl 修剪后残留的弱信号 |
| intracellular.ralstonia_like | 2,776 | 37 | 1.3% | GASRG/GTSRG 变异型 |
| intracellular.bacillus_like | 57 | 0 | 0% | 修剪后完全丢失 |
| extracellular.general | 509 | 0 | 0% | 修剪后完全丢失 |
| extracellular.lemoignei_like | 412 | 3 | 0.7% | 修剪后残留 |

**⚠ 重要说明**:
- 以上数据基于 **trimAl 修剪后**的比对，lipase box 是可变区域，已被 trimAl 大量去除
- 论文中 Lipase box 统计应以 **原始序列**（验证集 FASTA）的 motif 扫描结果为准:
  - intracellular: 10.1%
  - ralstonia: 9.8%
  - extracellular: 28.6%
  - extracellular_lemoignei: 53.2%
  - bacillus_type: 50.0%
- 这些百分比来自对完整序列的 G-X-S-X-G 正则匹配，不受 trimAl 影响

---

## 4. #16: Pfam 结构域注释

| 数据 | 说明 |
|---|---|
| 输入 | 6,532 条验证 PhaZ 蛋白序列 |
| 数据库 | Pfam-A (完整 release, ~18,000 HMM) |
| 方法 | HMMER hmmscan, E≤1e-5, 16 线程 |
| 策略 | **全部扫描**，无抽样 |

---

## 5. T4: 催化域保守性

| 数据 | 说明 |
|---|---|
| 输入 | c95 级别 trimAl 修剪后比对 (8,178 条) |
| 方法 | 位点保守性计算（香农熵 + consensus%） |
| 策略 | 全部使用 c95 比对，无需额外抽样 |

---

## 6. T5: 信号肽预测

| 数据 | 说明 |
|---|---|
| 输入 | 6,532 条验证 PhaZ 蛋白序列 |
| 方法 | Kyte-Doolittle 疏水性扫描（N端 30aa） |
| 策略 | **全部扫描**，无抽样 |

---

## 7. T6: 宏基因组种子序列

| 亚型 | 输入 | CD-HIT c70 | 最终种子 |
|---|---|---|---|
| intracellular.cupriavidus_like | 3,668 | CD-HIT 失败* | 10 |
| intracellular.ralstonia_like | 2,221 | CD-HIT 失败* | 10 |
| intracellular.bacillus_like | 5 | 全部保留 | 5 |
| extracellular.general | 404 | CD-HIT 失败* | 10 |
| extracellular.lemoignei_like | 234 | CD-HIT 失败* | 10 |
| **合计** | 6,532 | | **45** |

*CD-HIT 未在 PATH 中找到，回退方案：每亚型选最长 10 条序列

---

## 总结

| 分析 | 数据级别 | 抽样 |
|---|---|---|
| 核心统计（拷贝数/分布/亚型） | 验证集 6,532 | **无抽样** |
| 系统发育树 | c95 级别 8,178 | 各亚型独立建树（无抽样） |
| 联合树 | c95 抽稀 857 | 小亚型全保留，大亚型 ≤250 |
| 3D 结构 | 每亚型 5 条代表 | 最长序列 |
| 保守性 | c95 修剪后 | 无抽样 |
| 信号肽 | 验证集 | 无抽样 |
| Pfam 注释 | 验证集 | 无抽样 |
| 宏基因组种子 | 验证集 | 每亚型 10 条 |
