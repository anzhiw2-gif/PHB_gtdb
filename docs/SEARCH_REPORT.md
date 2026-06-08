# PHB_gtdb — PhaZ 搜索结果报告

**生成日期**: 2026-06-08  
**GTDB 版本**: Release R232  
**搜索方法**: Pyrodigal (in-process) + DIAMOND blastp  
**参考序列**: 14 条 NCBI 验证的 PhaZ/PHA depolymerase

---

## 1. 搜索概览

| 指标 | 古菌 (Archaea) | 细菌 (Bacteria) | 合计 |
|------|:---:|:---:|:---:|
| 搜索基因组数 | 10,122 | 165,468 | 175,590 |
| 含 PhaZ 基因组 | **2** | **7,068** | 7,070 |
| PhaZ 基因总数 | **2** | **16,486** | 16,488 |
| 命中率 | 0.02% | 4.27% | 4.03% |
| 覆盖门数 | 2 | 15 | 15 |
| 耗时 | 2h11m | 40h34m | 42h45m |

---

## 2. 古菌搜索结果

### ⚠️ 发现 2 个古菌含 PhaZ 同源序列

| 基因组 | 物种 | 古菌门 | 古菌纲 | Identity | E-value | 匹配参考 |
|---|---|---|---|---|---|---|
| GCA_003551945.1 | *Halovenus* sp. | Halobacteriota | Halobacteria | 71.4% | 4.71e-11 | WP_128854079.1 (Bacillus 胞内型) |
| GCA_043821585.1 | UBA73 sp. | Nanobdellota | Pacearchaeales | 52.5% | 5.45e-11 | WP_207907290.1 (P. lemoignei 胞外型) |

### 分析

- **Halovenus sp.** (盐古菌): 71.4% identity 匹配 Bacillus 胞内型 PhaZ — 高相似度提示可能通过**水平基因转移 (HGT)** 从细菌获得
- **UBA73 sp.** (DPANN 超门): 52.5% identity — 边缘 e-value (5.45e-11)，可能为保守 α/β 水解酶结构域的交叉匹配
- 古菌命中率 (0.02%) 极低，**总体上支持"古菌不携带 phaZ 基因"的结论**
- 建议: 对这两条序列做 reciprocal BLASTP + 结构域分析确认

---

## 3. 细菌搜索结果

### 3.1 门级分布

| 门 (Phylum) | 含 PhaZ 基因组数 | 占比 |
|---|---|---|
| **Pseudomonadota** | 6,681 | 94.5% |
| Actinomycetota | 244 | 3.5% |
| Bacillota | 51 | 0.7% |
| Myxococcota | 24 | 0.3% |
| Acidobacteriota | 11 | 0.2% |
| Bdellovibrionota | 10 | 0.1% |
| Bacteroidota | 10 | 0.1% |
| Chloroflexota | 9 | 0.1% |
| Desulfobacterota | 8 | 0.1% |
| Spirochaetota | 5 | 0.1% |
| Planctomycetota | 5 | 0.1% |
| Gemmatimonadota | 5 | 0.1% |
| Verrucomicrobiota | 3 | <0.1% |
| Thermotogota | 1 | <0.1% |
| Cyanobacteriota | 1 | <0.1% |

![门级分布](results/figures/phylum_distribution.png)

### 3.2 PhaZ 富集最多的基因组 (Top 10)

| 基因组 | PhaZ 数 | 最佳 Identity | 物种 | 分类 |
|---|---|---|---|---|
| GCF_017814975.1 | **12** | 67.6% | — | — |
| GCF_052113425.1 | **11** | **88.0%** | — | — |
| GCF_046535245.1 | **11** | 69.2% | — | — |
| GCF_017744395.1 | 11 | 43.3% | — | — |
| GCF_015767275.1 | 11 | 50.0% | — | — |
| GCA_963724595.1 | 11 | 50.0% | *Hyphomicrobium_A* | α-变形菌 |
| GCA_016177875.1 | 10 | **85.7%** | ***Paucimonas*** | β-变形菌 (参考来源) |
| GCA_050987995.1 | 10 | 50.5% | ***Cupriavidus pickettii*** | β-变形菌 (参考来源) |
| GCF_004768545.1 | 10 | 46.7% | — | — |
| GCF_046156345.1 | 10 | 43.9% | — | — |

> *斜体* = 参考序列来源属

### 3.3 参考序列匹配分布

| 参考序列 | 类型 | 匹配次数 | 代表物种 |
|---|---|---|---|
| BAA33394.1 | 胞内型 PhaZ1 | ★ 最多 | *C. necator* H16 |
| CAJ93939.1 | 胞内型 PhaZ2 | ★ 多 | *C. necator* H16 |
| CAJ95805.1 | 胞内型 PhaZ5 | ★ 多 | *C. necator* H16 |
| WKZ88401.1 | PHA depolymerase | ★ 多 | *R. pickettii* |
| UCA14981.1 | PHA depolymerase | ★ 多 | *R. pickettii* |
| BAA19791.1 | 胞外型 | 中等 | *D. acidovorans* |
| BAA35137.1 | 胞外型 | 中等 | *Acidovorax* sp. |
| AAA87070.1 | 胞外型 | 中等 | *C. testosteroni* |

> 胞内型 *C. necator* PhaZ 家族是最普遍的 PhaZ 类型，分布在多个细菌门中

---

## 4. 方法

### 基因预测

使用 [Pyrodigal](https://github.com/althonos/pyrodigal) v3.7.1 进行进程内基因预测（无 subprocess 开销）。

### 同源搜索

使用 [DIAMOND](https://github.com/bbuchfink/diamond) v2.1.24 blastp:

```
diamond blastp -q proteins.faa -d phaz_db.dmnd \
    -e 1e-10 --id 30 --query-cover 50 --max-target-seqs 5
```

### 筛选标准

| 参数 | 阈值 |
|---|---|
| E-value | ≤ 1e-10 |
| Identity | ≥ 30% |
| Query coverage | ≥ 50% |

### 14 条参考序列

详见 [`docs/PHAZ_REFERENCES.md`](docs/PHAZ_REFERENCES.md)，覆盖 8 属、3 门，包含胞外型 (8条) 和胞内型 (6条)。

---

## 5. 已知限制

1. **假单胞菌门偏倚**: 11/14 参考序列来自 Pseudomonadota，导致该门占命中 94.5%
2. **严格阈值**: id≥30% + qcov≥50% 可能遗漏远缘 PhaZ 同源物
3. **参考序列数有限**: 14 条无法覆盖 PhaZ 完整序列空间
4. **古菌命中待验证**: 2 个古菌命中需 reciprocal BLASTP + 结构域分析

---

## 6. 后续步骤

- [x] Step 1: PhaZ 基因搜索 ✅
- [ ] Step 2: 提取 PhaZ 蛋白序列
- [ ] Step 3: MAFFT 多序列比对 + trimAl
- [ ] Step 4: IQ-TREE 系统发育分析
- [ ] Step 5: 蛋白结构域注释 (Pfam/InterPro)
- [ ] Step 6: 可视化 (门分布、系统发育树、热图)

---

**报告生成时间**: 2026-06-08 18:00 CST  
**分析服务器**: T141 (10.16.1.141)  
**GitHub**: [anzhiw2-gif/PHB_gtdb](https://github.com/anzhiw2-gif/PHB_gtdb)
