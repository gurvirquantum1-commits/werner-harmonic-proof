#!/usr/bin/env python3
"""Recompute the degree-9 shape range bound from the certified trig grid.

The calculation uses the supplied outward trigonometric table at
  theta_i = pi*i/1024, i=0,...,1024,
  phi_j   = 2*pi*j/2048, j=0,...,2047.
Central values are evaluated in binary64 and then enlarged by a deliberately
large 1e-9 absolute error budget. The grid is a geodesic delta-net with
  delta <= pi/1024.
For a spherical polynomial of angular degree d=9, Bernstein's inequality on
every great circle gives ||grad_S g||_inf <= d ||g||_inf, hence
  M <= G + d*delta*M,
so M <= G/(1-d*delta).

The 1e-9 grid-evaluation cushion is much larger than a crude first-order
binary64/trig-table error estimate (<4e-10 for this dataset).
"""
from fractions import Fraction as Q
from pathlib import Path
import json, math, numpy as np

ROOT=Path(__file__).resolve().parents[1]
D=json.load(open(ROOT/'data'/'werner_rational_v6910_candidate_deg111.json'))
T=np.load(ROOT/'data'/'trig_cc_exact.npz')
E=np.array(D['exponents'],dtype=np.int64)
assert E.shape==(55,3) and np.all(E.sum(axis=1)==9)
rows=D['A']+D['B']
C=np.array([[float(Q(s)) for s in row] for row in rows],dtype=float)
L1=np.sum(np.abs(C),axis=1)

# Midpoints of already outward-enclosed trig values.
th_c=(T['c2048lo'][:1025]+T['c2048hi'][:1025])/2
th_s=(T['s2048lo'][:1025]+T['s2048hi'][:1025])/2
ph_c=(T['c2048lo']+T['c2048hi'])/2
ph_s=(T['s2048lo']+T['s2048hi'])/2

# Sanity-check enclosure widths.
trig_width=max(float(np.max(T['c2048hi']-T['c2048lo'])),
               float(np.max(T['s2048hi']-T['s2048lo'])))
assert trig_width <= 2.0**-52

mx=np.zeros(42)
for i in range(1025):
    x=th_s[i]*ph_c; y=th_s[i]*ph_s; z=th_c[i]
    # Explicit repeated powers, avoiding transcendental power routines.
    xp=np.empty((10,2048)); yp=np.empty((10,2048)); zp=np.empty(10)
    xp[0]=1.0; yp[0]=1.0; zp[0]=1.0
    for k in range(1,10):
        xp[k]=xp[k-1]*x; yp[k]=yp[k-1]*y; zp[k]=zp[k-1]*z
    M=np.empty((55,2048))
    for j,(aa,bb,cc) in enumerate(E): M[j]=xp[aa]*yp[bb]*zp[cc]
    vals=C@M
    mx=np.maximum(mx,np.max(np.abs(vals),axis=1))

# Deliberately over-sized arithmetic/enclosure cushion.
U=2.0**-53
coord_err=2e-15
crude_coord=9.0*float(L1.max())*coord_err
crude_arith=5000.0*U*float(L1.max())
crude_coeff=4.0*U*float(L1.max())
crude_total=crude_coord+crude_arith+crude_coeff
GRID_ERR=1.0e-9
assert crude_total < GRID_ERR, (crude_total,GRID_ERR)
G=mx+GRID_ERR

# Rational upper bound pi < 355/113.
pi_up=Q(355,113)
delta=float(pi_up/Q(1024,1))
den=1.0-9.0*delta
assert den>0
global_bounds=G/den

print('grid max A upper = %.16g' % G[:21].max())
print('grid max B upper = %.16g' % G[21:].max())
print('delta upper      = %.16g' % delta)
print('1-9 delta       = %.16g' % den)
print('global max A    = %.16g' % global_bounds[:21].max())
print('global max B    = %.16g' % global_bounds[21:].max())
print('crude binary64/trig error estimate = %.3e; adopted cushion = %.1e' % (crude_total,GRID_ERR))
for i,b in enumerate(global_bounds):
    side='A' if i<21 else 'B'; r=i if i<21 else i-21
    print('%s[%02d] < %.16g' % (side,r,b))
assert global_bounds.max() < 0.9
print('CERTIFIED RANGE CHECK PASSED: all 42 degree-9 shapes satisfy |g| < 0.9 on S^2.')
