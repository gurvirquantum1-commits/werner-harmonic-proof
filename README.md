# Werner harmonic local model

Verification code and exact candidate data for *A harmonic approach to local models for two-qubit Werner states*, by Gurvir Singh.

Repository: [gurvirquantum1-commits/werner-harmonic-proof](https://github.com/gurvirquantum1-commits/werner-harmonic-proof).

The certificate supports the projective-measurement visibility bound

$$v_c^{\mathrm{proj}} \geq \frac{6914213000}{10002162779} = 0.69127179318824\ldots.$$

This is a lower bound, not the exact critical visibility.

The full run was reproduced on 29 October 2026. See the [fresh verification record](recorded_results/reproduction_2026-10-29/VERIFICATION.md) for results, environment and logs.

## Verify

Use Python 3.11 or later. From this directory:

```sh
python -m venv .venv
```

Activate the environment using the command for your system:

```powershell
# Windows PowerShell
.venv\Scripts\Activate.ps1
```

```sh
# macOS / Linux
source .venv/bin/activate
```

Then run:

```sh
python -m pip install -r requirements.txt
python run_audit.py
```

Package installation needs internet access; verification itself runs offline. Do not use `python -O` or `python -OO`: the checks use assertions. Allow several minutes or longer depending on hardware, and at least 2 GB of free disk space for regenerated arrays.

A successful run ends with:

```text
Manuscript exact shortened bound: 6914213000/10002162779
ALL FULL VERIFICATION CHECKS PASSED
```

The command checks file hashes, exact rational weights and response coefficients, the full spherical range grid and Bernstein extension, and stored-certificate arithmetic. It then independently rebuilds trigonometric tables, Legendre functions, harmonic projections and tail bounds, checks the archived enclosures, and runs the original regeneration route using its NumPy compatibility implementation. Numba is not required.

Fresh arrays and certificates go to `work/`; per-stage logs and environment/timing information go to `logs/`. These generated folders are excluded from Git. A failed stage stops the run and names its log.

## Files

| Path | Contents |
| --- | --- |
| `run_audit.py` | Complete verification entry point |
| `data/` | Supplied exact candidate, tables and certificate |
| `scripts/` | Exact checks and both numerical reconstruction routes |
| `recorded_results/` | Historical audit records and fresh reproduction record |
| `REPORT.md` | Supplied audit report, mathematical details and numerical error assumptions |
| `SHA256.json` | Integrity manifest for distributed files |

For short checks alone, run `python scripts/verify_all.py`. **That command does not reconstruct the harmonic certificate**; use `python run_audit.py` for full verification.

## Scope

The result concerns arbitrary single-copy projective measurements on two-qubit Werner states. Numerical enclosures assume standard IEEE-754 arithmetic; this is computer-assisted verification, not proof-assistant formalization. The final shortened fraction is checked with exact rational arithmetic. Auxiliary local correction kernels exist by Grothendieck's inequality; explicit deterministic strategy tables for them are not supplied.

The supplied `REPORT.md` and `data/certificate_summary.json` are retained as historical inputs. The summary's statement that the raw generator is absent is stale: its component scripts are included. The report describes the earlier seven-stage independent run; this repository's full command additionally runs the NumPy regeneration and its certificate check. Historical logs are reference material, not substitutes for recomputation.
