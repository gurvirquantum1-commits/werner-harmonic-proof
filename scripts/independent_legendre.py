"""High-precision normalized associated Legendre recurrence at exact CC nodes.
For each m,l, a global ball radius is propagated with a 1e-90 local allowance.
At precision110 and |intermediates|<1000 this exceeds Decimal arithmetic,
constant, node and normalization errors by >1e10. Bounded exact functions
satisfy |Pbar_lm|<=sqrt((2l+1)/2)<13. Each stored float is widened outward.
"""
from decimal import Decimal as D,getcontext
from pathlib import Path
import numpy as np,json,time
getcontext().prec=110
root=Path(__file__).resolve().parents[1];work=root/'work'
nodes=json.load(open(work/'independent_decimal_nodes.json'));z=np.array([D(x) for x in nodes['z']],dtype=object);st=np.array([D(x) for x in nodes['st']],dtype=object)
EPS=D('1e-90');L=151;N=len(z);diag=np.full(N,1/D(2).sqrt(),dtype=object);dr=EPS
P=np.memmap(work/'independent_Pbar.dat',mode='w+',dtype=np.float64,shape=(L+1,L+1,N));rad=np.zeros((L+1,L+1));mr=D(0);t=time.time()
for m in range(L+1):
 if m:
  fac=-(D(1)+D(1)/D(2*m)).sqrt();diag=fac*st*diag;dr=abs(fac)*dr+EPS
 p=diag.copy();prev=None;r=dr;rp=D(0)
 for l in range(m,L+1):
  p_float=np.array(p,dtype=np.float64);P[m,l]=p_float
  # Conversion to nearest binary64 gives absolute error <= 2e-15 for |P|<13.
  rad[m,l]=float(r)+2e-15
  mr=max(mr,r)
  assert r<D('1e-15') and np.max(np.abs(p_float))<13
  if l==L:break
  if l==m:
   A=D(2*m+3).sqrt();nxt=A*z*p;rn=A*r+EPS
  else:
   ll=l+1;A=(D(4*ll*ll-1)/D(ll*ll-m*m)).sqrt();B=(D((2*ll+1)*((ll-1)*(ll-1)-m*m))/D((2*ll-3)*(ll*ll-m*m))).sqrt()
   nxt=A*z*p-B*prev;rn=A*r+B*rp+EPS
  prev,p=p,nxt;rp,r=r,rn
 if m%20==0:print('m',m,'seconds',time.time()-t,'max_decimal_radius',mr,flush=True)
P.flush();np.save(work/'independent_Pbar_radius.npy',rad)
print('DONE maximum Decimal ball radius',mr,'seconds',time.time()-t,flush=True)
