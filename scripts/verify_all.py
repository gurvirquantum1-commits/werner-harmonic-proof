#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys
HERE=Path(__file__).resolve().parent
scripts=['verify_candidate_integrity.py','verify_response_exact.py','verify_range_grid_bernstein.py','verify_harmonic_certificate.py']
for s in scripts:
    print('\n'+'='*78+'\nRUN',s+'\n'+'='*78,flush=True)
    subprocess.run([sys.executable,str(HERE/s)],check=True)
print('\nALL CERTIFICATE CHECKS PASSED.')
