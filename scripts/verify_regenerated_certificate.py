#!/usr/bin/env python3
"""Check that the raw regeneration independently clears the manuscript theorem."""
import json, os
from fractions import Fraction as F
HERE=os.path.dirname(os.path.abspath(__file__))
ROOT=os.path.dirname(HERE); WORK=os.path.join(ROOT,'work')
r=json.load(open(os.path.join(WORK,'werner_harmonic_certificate_deg111_regenerated.json')))
# Deliberately shortened exact manuscript endpoints.
v0q=F(6914213,10_000_000); Dq=F(1213,10_000_000); K=F(1783,1000)
assert r['v0_lower'] > float(v0q)
assert r['D_total_upper'] < float(Dq)
exact=v0q/(1+K*Dq)
assert exact == F(6914213000,10002162779)
assert r['v_cert_lower'] > float(exact)
print('raw regenerated v0 lower:',r['v0_lower'])
print('raw regenerated D upper:',r['D_total_upper'])
print('raw regenerated corrected visibility:',r['v_cert_lower'])
print('manuscript exact shortened theorem:',exact,'=',float(exact))
print('PASS: independent low-level regeneration clears the manuscript theorem.')
