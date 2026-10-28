#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
D=json.load(open(ROOT/'data'/'werner_rational_v6910_candidate_deg111.json'))
E=[tuple(x) for x in D['exponents']]
assert len(E)==55 and len(set(E))==55 and all(sum(e)==9 for e in E)
assert len(D['A'])==len(D['B'])==len(D['weights'])==21
assert all(len(row)==55 for row in D['A']+D['B'])
# Every decimal/rational string parses exactly.
for row in D['A']+D['B']:
    [Q(s) for s in row]
w=[Q(s) for s in D['weights']]
assert all(x>0 for x in w)
assert sum(w,Q(0))==1
assert D['response']['degree']==111
assert len(D['response']['cheb_coeffs'])==112
print('CANDIDATE INTEGRITY PASSED')
print('21 positive exact rational weights sum exactly to 1.')
print('42 homogeneous degree-9 rational shapes, 55 monomials each.')
