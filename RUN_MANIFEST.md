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
