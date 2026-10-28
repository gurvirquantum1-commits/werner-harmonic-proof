#!/usr/bin/env python3
"""Generate outward-enclosed trigonometric values and Clenshaw--Curtis weights.

This is the first stage of the independent low-level regeneration of the
ell<=151 harmonic certificate.  The angular table is obtained by rigorous half-angle refinement of the
supplied outward 2048-angle trigonometric table.  The CC weights use the
standard even-N cosine formula, with outward trig intervals and an explicit
floating-point summation envelope.
"""
import os, math
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__))
ROOT=os.path.dirname(HERE)
WORK=os.path.join(ROOT,'work')
os.makedirs(WORK,exist_ok=True)
INF=np.inf
U=2.0**-53

def down(x): return np.nextafter(x,-INF)
def up(x): return np.nextafter(x,INF)

# Rigorous half-angle refinement of the supplied certified 2048-angle table.
# data/trig_cc_exact.npz contains outward enclosures of
# cos(2*pi*j/2048), sin(2*pi*j/2048), and pi.  For 0<=j<=2048 we need
# theta_j=2*pi*j/4096.  Since cos(2 theta_j) is already certified, use
#   |cos theta|=sqrt((1+cos 2theta)/2),
#   sin theta =sqrt((1-cos 2theta)/2)  (theta in [0,pi]),
# with the known sign of cos theta.  Directed nextafter after sqrt keeps
# the binary64 intervals outward.  The second half follows by symmetry.
base=np.load(os.path.join(ROOT,'data','trig_cc_exact.npz'))
c2lo=base['c2048lo']; c2hi=base['c2048hi']
N4=4096
cl=np.empty(N4+1);ch=np.empty(N4+1);sl=np.empty(N4+1);sh=np.empty(N4+1)
def sqrt_iv(lo,hi):
    lo=max(0.0,float(lo)); hi=max(0.0,float(hi))
    return down(math.sqrt(lo)), up(math.sqrt(hi))
for j in range(N4//2+1):
    if j==2048:
        c2l=c2h=1.0  # cos(2*pi)
    else:
        c2l=float(c2lo[j]);c2h=float(c2hi[j])
    aclo,achi=sqrt_iv((1+c2l)/2,(1+c2h)/2)
    sloj,shij=sqrt_iv((1-c2h)/2,(1-c2l)/2)
    if j<=1024:
        cl[j],ch[j]=aclo,achi
    else:
        cl[j],ch[j]=-achi,-aclo
    sl[j],sh[j]=sloj,shij
for j in range(N4//2+1,N4+1):
    k=N4-j
    cl[j],ch[j]=cl[k],ch[k]
    sl[j],sh[j]=-sh[k],-sl[k]
for j,c,sn in [(0,1.,0.),(1024,0.,1.),(2048,-1.,0.),(3072,0.,-1.),(4096,1.,0.)]:
    cl[j]=ch[j]=c;sl[j]=sh[j]=sn
np.savez(os.path.join(WORK,'trig4096_mp_iv.npz'),cl=cl,ch=ch,sl=sl,sh=sh)

# Clenshaw--Curtis weights for N=2048.  For 1<=j<N,
# w_j=(2/N)[1-sum_{k=1}^{N/2-1} 2 cos(2j k pi/N)/(4k^2-1)
#                - (-1)^j/(N^2-1)].
N=2048; NZ=N+1
js=np.arange(1,N,dtype=np.int64); ks=np.arange(1,N//2,dtype=np.int64)
idx=(2*(js[:,None]*ks[None,:]))%N4
co_lo=cl[idx];co_hi=ch[idx]
a=np.array([2.0/(4*int(k)*int(k)-1) for k in ks])
TL=down(co_lo*a[None,:]); TH=up(co_hi*a[None,:])
sl0=np.sum(TL,axis=1); sh0=np.sum(TH,axis=1)
gam=(2*len(ks)*U)/(1-2*len(ks)*U)
sl0=down(sl0-up(gam*np.sum(np.abs(TL),axis=1)))
sh0=up(sh0+up(gam*np.sum(np.abs(TH),axis=1)))
end_lo=down(1.0/(N*N-1));end_hi=up(1.0/(N*N-1))
last_lo=np.where(js%2==0,end_lo,-end_hi);last_hi=np.where(js%2==0,end_hi,-end_lo)
br_lo=down(1.0-sh0-last_hi); br_hi=up(1.0-sl0-last_lo)
fac_lo=down(2.0/N);fac_hi=up(2.0/N)
wlo=np.empty(NZ);whi=np.empty(NZ)
wlo[0]=wlo[-1]=end_lo;whi[0]=whi[-1]=end_hi
wlo[1:N]=down(br_lo*fac_lo);whi[1:N]=up(br_hi*fac_hi)
np.savez(os.path.join(WORK,'cc2048_weights_iv.npz'),lo=wlo,hi=whi)
print('trig max width:',max(np.max(ch-cl),np.max(sh-sl)))
print('CC weight sum enclosure:',down(np.sum(wlo)),up(np.sum(whi)))
print('CC max weight width:',np.max(whi-wlo))
