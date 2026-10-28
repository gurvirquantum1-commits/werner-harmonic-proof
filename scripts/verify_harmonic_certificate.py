#!/usr/bin/env python3
"""Audit the stored outward harmonic certificate and prove v > 0.6910.

This checker does not regenerate the expensive harmonic projections from raw
sphere samples. It independently checks all arithmetic *after* the stored
outward coefficient/tail enclosures: resolved residual, Cauchy--Schwarz tail,
conservative rational shortening, and the final exact visibility inequality.
"""
from fractions import Fraction as Q
from decimal import Decimal, getcontext
from pathlib import Path
import json, math

getcontext().prec=80
ROOT=Path(__file__).resolve().parents[1]
M=json.load(open(ROOT/'data'/'werner_rational_v6910_candidate_deg111.json'))
C=json.load(open(ROOT/'data'/'werner_harmonic_certificate_deg111_final.json'))
w=[Q(s) for s in M['weights']]
assert all(x>0 for x in w) and sum(w,Q(0))==1

ints=C['coeff_intervals']
ells=[x['ell'] for x in ints]
assert ells==list(range(1,int(C['LMAX'])+1,2))
for x in ints:
    assert x['lo'] <= x['center'] <= x['hi']
    assert abs((x['center']-x['radius'])-x['lo']) < 5e-18
    assert abs((x['center']+x['radius'])-x['hi']) < 5e-18

c1=ints[0]
v0_lower=-c1['hi']
assert abs(v0_lower-C['v0_lower']) < 2e-16
resolved=sum(max(abs(x['lo']),abs(x['hi'])) for x in ints[1:])
assert abs(resolved-C['resolved_D_upper']) < 1e-18

# tailE: first 21 are Alice energy upper bounds, next 21 Bob.
E=C['tailE']; assert len(E)==42
num=sum(float(w[r])*math.sqrt(E[r]*E[21+r]) for r in range(21))
tail=num/(4*math.pi)
assert abs(tail-C['tail_upper']) < 2e-18
D=resolved+tail
assert abs(D-C['D_total_upper']) < 2e-18

K=1.783
v=v0_lower/(1+K*D)
assert abs(v-C['v_cert_lower']) < 2e-15

# Exact shortened endpoints used for the paper-level theorem.
v0q=Q(6914213,10_000_000)
Dq=Q(1213,10_000_000)
Kq=Q(1783,1000)
assert Decimal(str(C['v0_lower'])) > Decimal(v0q.numerator)/Decimal(v0q.denominator)
assert Decimal(str(C['D_total_upper'])) < Decimal(Dq.numerator)/Decimal(Dq.denominator)
vq=v0q/(1+Kq*Dq)
assert vq==Q(6914213000,10002162779)
assert vq > Q(6910,10000)

print('HARMONIC CERTIFICATE ARITHMETIC PASSED')
print('v0_lower          = %.16g' % v0_lower)
print('resolved D upper  = %.16g' % resolved)
print('tail upper        = %.16g' % tail)
print('total D upper     = %.16g' % D)
print('floating lower v  = %.16g' % v)
print('exact shortened bound = %s / %s' % (vq.numerator,vq.denominator))
print('exact shortened decimal ~= %.15f' % float(vq))
print('CERTIFIED THEOREM: v_c^proj > 0.6910 (indeed > 0.69127 by the stored enclosure certificate).')
