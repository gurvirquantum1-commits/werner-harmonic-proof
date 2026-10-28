"""Direct polynomial evaluation + direct DFT, independently of supplied FFT
and 32-point shape reconstruction. Uses independently regenerated tables.
Error derivation: coordinates within 5u, degree9 monomial absolute error <60u,
coefficient rounding <u|c| and length55 dot <56u*sum|c|.
200u*max row L1 covers their sum with >60% reserve. The odd Chebyshev
recurrence has |U_k(x)|<=k+1; 40u per recurrence forcing gives weighted
roundoff <1e-12 for this coefficient list. We adopt 2e-11.
Direct DFT dot length2048 roundoff <gamma_4096 * 2pi <3e-12;
we adopt 1e-10 (including phase/scale and coefficient errors).
"""
from pathlib import Path
from fractions import Fraction as F
import numpy as np,json,math,time
root=Path(__file__).resolve().parents[1];work=root/'work'
D=json.load(open(root/'data/werner_rational_v6910_candidate_deg111.json'));E=np.array(D['exponents']);C=np.array([[float(F(s)) for s in row] for row in D['A']+D['B']]);co=np.array([float(F(s)) for s in D['response']['cheb_coeffs']])[1::2]
T=np.load(work/'independent_trig4096.npz');W=np.load(work/'independent_cc2048.npz');nodes=json.load(open(work/'independent_decimal_nodes.json'));pi=float(nodes['pi']);dp=2*pi/2048
z=(T['cl'][:2049]+T['ch'][:2049])/2;st=(T['sl'][:2049]+T['sh'][:2049])/2
cm=(T['cl'][::2][:2048]+T['ch'][::2][:2048])/2;sm=(T['sl'][::2][:2048]+T['sh'][::2][:2048])/2
j=np.arange(2048);m=np.arange(152);ix=(j[:,None]*m[None,:])%2048;TC=cm[ix];TS=-sm[ix]
Fout=np.memmap(work/'independent_Fmid.dat',mode='w+',dtype=np.complex128,shape=(152,2049,42));ener=np.zeros((2049,42));maxg=0.;maxq=0.;t=time.time()
for i in range(2049):
 x=st[i]*cm;y=st[i]*sm
 xp=np.ones((10,2048));yp=np.ones((10,2048));zp=np.ones(10)
 for k in range(1,10):xp[k]=xp[k-1]*x;yp[k]=yp[k-1]*y;zp[k]=zp[k-1]*z[i]
 mono=np.array([xp[a]*yp[b]*zp[c] for a,b,c in E]);g=C@mono;maxg=max(maxg,float(abs(g).max()))
 u=g/0.9;xx=2*u*u-1;r0=u;q=co[0]*r0;r1=u*(2*xx-1);q+=co[1]*r1
 for k in range(2,len(co)):
  r2=2*xx*r1-r0;q+=co[k]*r2;r0,r1=r1,r2
 maxq=max(maxq,float(abs(q).max()));Fout[:,i,:]=(dp*(q@TC+1j*(q@TS))).T
 ener[i]=dp*np.sum(q*q,axis=1)
 if i%256==0:print(i,'seconds',time.time()-t,flush=True)
Fout.flush();total=((W['lo']+W['hi'])/2)@ener
u=2**-53;gerr=200*u*float(abs(C).sum(axis=1).max());lip=float(sum(abs(F(s))*k*k for k,s in enumerate(D['response']['cheb_coeffs'])))/.9;eta_derived=lip*gerr+2e-11
assert eta_derived<2e-8
np.savez(work/'independent_meta.npz',total=total,gerr=gerr,lip=lip,eta_derived=eta_derived,eta_f=2e-8,max_g=maxg,max_q=maxq)
print('DONE',time.time()-t,'gerr',gerr,'eta_derived',eta_derived,'maxg',maxg,'maxq',maxq,flush=True)
