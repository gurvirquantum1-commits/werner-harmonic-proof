#!/usr/bin/env python3
"""Exact rational verification of the degree-111 soft-sign response.

This script reconstructs F_m(a x) exactly for m=550, a=9/10,
converts it to the Chebyshev basis, recomputes the omitted tail after
degree 111, and checks that the stored response coefficients are exactly
the scaled truncation T/(1+eps).

No floating-point arithmetic is used for the identities.
"""
from fractions import Fraction as Q
from pathlib import Path
import json, math

ROOT = Path(__file__).resolve().parents[1]
D = json.load(open(ROOT/'data'/'werner_rational_v6910_candidate_deg111.json'))
R = D['response']
m = int(R['m_source'])
a = Q(R['domain_bound'])
deg = int(R['degree'])

# F_m(x) = [int_0^x (1-t^2)^m dt] / [int_0^1 (1-t^2)^m dt].
den = sum(((-1)**j) * Q(math.comb(m,j), 2*j+1) for j in range(m+1))
assert den > 0

# Exact Chebyshev coefficients of F_m(a y).
# For odd n: y^n = 2^(1-n) sum_{r=0}^{(n-1)/2} C(n,r) T_{n-2r}(y).
cheb = [Q(0) for _ in range(2*m+2)]
for j in range(m+1):
    n = 2*j + 1
    mono = ((-1)**j) * Q(math.comb(m,j), 2*j+1) * (a**n) / den
    fac = Q(1, 2**(n-1))
    for r in range((n-1)//2 + 1):
        k = n - 2*r
        cheb[k] += mono * fac * math.comb(n,r)

# Uniform truncation bound because |T_k(y)| <= 1 for |y|<=1.
eps = sum(abs(cheb[k]) for k in range(deg+1, len(cheb)))
stored_eps = Q(R['tail_epsilon_exact'])
assert eps == stored_eps, 'stored tail epsilon does not match exact reconstruction'
assert abs(float(eps)-float(R['tail_epsilon_decimal'])) < 1e-18

scale = Q(R['scale_exact'])
assert scale == Q(1,1)/(Q(1,1)+eps), 'scale must be exactly 1/(1+epsilon)'

stored = [Q(s) for s in R['cheb_coeffs']]
assert len(stored) == deg+1
for k in range(deg+1):
    assert stored[k] == cheb[k]*scale, f'Chebyshev coefficient mismatch at k={k}'
    if k % 2 == 0:
        assert stored[k] == 0, f'even coefficient {k} should vanish'

print('EXACT RESPONSE CERTIFICATE PASSED')
print(f'm = {m}, domain |g| <= {a}, retained Chebyshev degree = {deg}')
print('epsilon_exact =', eps)
print('epsilon_decimal ~= %.17g' % float(eps))
print('scale_exact = 1/(1+epsilon) verified')
print('all stored Chebyshev coefficients match exact scaled truncation')
print('Therefore, if |g|<=9/10 then |f(g)|<=1 and q=(1+f)/2 lies in [0,1].')
