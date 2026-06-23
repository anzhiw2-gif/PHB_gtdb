# PHB_gtdb Run Manifest

## Runtime Location

- Server: T141 (`10.16.1.141`)
- Hostname observed: `root123-PowerEdge-T630`
- Project path: `/home/data/haoyu/PHB_gtdb`
- GTDB root: `/home/data/haoyu/GTDB`
- GTDB release: R232

## Expected Core Result Counts

| File | Expected count |
|---|---:|
| `data/processed/phaz_proteins_all.fasta` | 8,768 |
| `data/processed/phaz_proteins_c95.fasta` | 8,731 |
| `data/processed/phaz_proteins_filtered.fasta` | 7,478 |
| `data/processed/phaz_proteins_validated.fasta` | 6,532 |
| `data/processed/phaz_bacillus_type.fasta` | 60 |
| `data/processed/phaz_extracellular.fasta` | 596 |
| `data/processed/phaz_extracellular_lemoignei.fasta` | 434 |
| `data/processed/phaz_intracellular.fasta` | 4,693 |
| `data/processed/phaz_ralstonia.fasta` | 2,948 |
| `data/processed/phaz_bacillus_type_trim.fasta` | 57 |
| `data/processed/phaz_extracellular_trim.fasta` | 509 |
| `data/processed/phaz_extracellular_lemoignei_trim.fasta` | 412 |
| `data/processed/phaz_intracellular_trim.fasta` | 4,424 |
| `data/processed/phaz_ralstonia_trim.fasta` | 2,776 |

| Search result | Expected rows | Expected `phaZ_count` sum |
|---|---:|---:|
| `data/processed/phb_search_results.tsv` | 7,068 | 16,486 |
| `data/processed/archaea_phb_search_results.tsv` | 2 | 2 |

## Validation

Run the standard result check from the project root:

```bash
python3 scripts/check_results.py
```

The checker verifies FASTA counts, search TSV counts, tree files, and figure files.

## Notes

- The server runtime directory may not itself be a git repository. Keep this manifest together with the scripts used to generate results.
- `phaz_proteins_c95.fasta` is the canonical c95 deduplicated file. Avoid using the ambiguous `phaz_proteins_dedup.fasta` name for final reporting.
- Historical failed logs may exist in `results/logs/`; use the latest successful logs and `scripts/check_results.py` for final validation.

## Tracked Figure Archive

The repository tracks the lightweight source data and final Nature-style figures:

| Directory | Purpose |
|---|---|
| `figure_data/` | Lightweight source data used by `scripts/07_nature_figures.py` |
| `figures/nature/` | Final Figure 1-5 outputs in PDF, PNG, and SVG |
| `figures/nature/source_data/` | Source data copied with the generated figures |

Expected figure files:

```text
figure1_workflow_funnel.pdf/png/svg
figure2_phylum_heatmap.pdf/png/svg
figure3_subtype_lipase.pdf/png/svg
figure4_genera_phylogeny.pdf/png/svg
figure5_grodon_growth_comparison_hmm_allmatched.pdf/png/svg
```

Figure titles and captions are maintained in `docs/FIGURE_CAPTIONS.md`.

## gRodon2 Growth-Rate Extension

Expert-suggested extension:

```text
Pyrodigal CDS/protein prediction
  -> Pfam ribosomal protein HMM + HMMER
  -> ribosomal proteins as highly expressed genes
  -> gRodon2 predictGrowth
  -> same-genus phaZ+ vs phaZ- comparison
```

Key files:

| File | Purpose |
|---|---|
| `scripts/08_prepare_ribosomal_hmms.py` | Build a reproducible Pfam ribosomal protein HMM subset |
| `scripts/08_grodon_growth.py` | Build matched manifest, predict CDS/HEG, call gRodon2, summarize same-genus differences |
| `scripts/08_run_grodon_one.R` | Single-genome `gRodon::predictGrowth` wrapper |
| `scripts/09_monitor_grodon_progress.py` | Live terminal/PNG monitor for the all-matched gRodon2 run |
| `scripts/10_balance_grodon_by_genus.py` | Create strict same-genus balanced `phaZ+`/`phaZ-` tables and failed-genome summaries |
| `scripts/11_grodon_growth_stats.py` | Run genus-level statistical tests and generate Figure 5 |
| `docs/GROWTH_RATE_ANALYSIS.md` | Design, pilot results, commands, and interpretation boundaries |
| `docs/GRODON2_FINAL_STATS.md` | Final balanced statistics, Figure 5 interpretation, and result boundaries |

Pilot status:

```text
Strict Pfam HMM pilot: 12/12 genomes succeeded
n_highly_expressed: 35-45
```

Formal all-matched manifest:

```text
8,788 genomes
4,394 phaZ+
4,394 same-genus phaZ-
903 genera
```

Formal all-matched run completed on T141:

```text
8,788 genomes processed
8,735 status=ok
53 failed
3 status=ok records lacked a valid growth-rate value
8,692 genomes retained after strict genus-level balancing
4,346 phaZ+
4,346 phaZ-
899 genera retained for genus-level tests
```

Final statistical result:

```text
Wilcoxon signed-rank P = 0.459
Exact sign test P = 0.386
Within-genus permutation P = 0.670
Conclusion: no significant global growth-rate difference between phaZ+ and phaZ- genomes after strict same-genus balancing.
```
