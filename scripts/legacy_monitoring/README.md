# Legacy Monitoring Helpers

These shell scripts were used during the early GTDB download, extraction, and
Step 1 search phases. They are kept for provenance but are not required for the
current reproducible workflow.

Current monitoring entry point:

```bash
python scripts/09_monitor_grodon_progress.py --once --no-clear
```

Legacy scripts:

| Script | Historical purpose |
|---|---|
| `download_monitor.sh` | Monitor GTDB R232 archive download |
| `extract_monitor.sh` | Monitor GTDB R232 archive extraction |
| `monitor_step1_search.sh` | Monitor the original bacterial Step 1 search |
| `monitor_search.sh` | Monitor bacterial/archaeal search progress |

