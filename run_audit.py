#!/usr/bin/env python3
"""Reproduce the independent Werner certificate audit using Python + NumPy."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,platform
ROOT=Path(__file__).resolve().parent
if sys.flags.optimize:raise RuntimeError('Run without -O or -OO: verification uses assertions.')
import numpy as np
manifest=json.load(open(ROOT/'SHA256.json'))
for name,want in manifest.items():
 got=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
 if got!=want:raise RuntimeError('Packaged file checksum mismatch: '+name)
(ROOT/'work').mkdir(exist_ok=True);(ROOT/'logs').mkdir(exist_ok=True)
env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);env['OPENBLAS_NUM_THREADS']='2';env['OMP_NUM_THREADS']='2';env['PYTHONUTF8']='1'
steps=['verify_all.py','generate_quadrature_tables.py','independent_tables.py','independent_legendre.py','independent_direct_compute.py','independent_certify.py','independent_crosschecks.py','compute_harmonic_data_numpy.py','certify_harmonics_raw.py','verify_regenerated_certificate.py']
summary={'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'stages':[],'passed':False};started=time.perf_counter()
for script in steps:
 print('RUN',script,flush=True)
 before=time.perf_counter()
 with open(ROOT/'logs'/(script+'.log'),'w',encoding='utf-8') as log:
  r=subprocess.run([sys.executable,'-u',str(ROOT/'scripts'/script)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=env)
 summary['stages'].append({'script':script,'returncode':r.returncode,'seconds':round(time.perf_counter()-before,3)})
 (ROOT/'logs/run_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
 if r.returncode:raise RuntimeError('Failed '+script+'; see logs/'+script+'.log')
 print('PASS',script,flush=True)
cert=json.load(open(ROOT/'work/independent_certificate.json'))
print('Independent conservative visibility lower bound:',cert['v_cert_lower'])
summary.update(passed=True,seconds=round(time.perf_counter()-started,3),independent_visibility_lower=cert['v_cert_lower'],exact_shortened_bound='6914213000/10002162779')
(ROOT/'logs/run_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print('Manuscript exact shortened bound: 6914213000/10002162779')
print('ALL FULL VERIFICATION CHECKS PASSED')
