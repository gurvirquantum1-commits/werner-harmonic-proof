# Reproduction — 29 October 2026

All ten stages of `python run_audit.py` passed on Windows with Python 3.12.14 and NumPy 2.3.5. Total elapsed time: 274.5 seconds.

| Route | Visibility lower bound |
| --- | ---: |
| Independent direct reconstruction | 0.691275623770854 |
| Original route, NumPy implementation | 0.691275621299138 |
| Manuscript shortened fraction | 6914213000/10002162779 |

The independent run checked containment of all 76 archived harmonic intervals and all 42 archived tail-energy bounds. Exact candidate, scalar-response, range and final rational checks passed. The two routes use different numerical implementations and error budgets, so their bounds differ slightly.

This records execution of the supplied checks and reconstruction algorithms under their documented numerical assumptions. It is not a proof-assistant formalization.


Fresh per-stage logs, generated certificate JSON files and `run_summary.json` are stored beside this record. Earlier results elsewhere in `recorded_results/` are historical supplied records.
