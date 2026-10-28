"""Independent Decimal/Taylor regeneration; no uploaded trig values used as inputs.
All Decimal work at 110 significant digits. Machin arctan and sin/cos series
are truncated only after next term <1e-115. We use absolute 1e-95 radius,
which dominates Taylor remainders and <1e6 Decimal roundings at this precision.
Binary64 output endpoints are widened one nextafter step after conversion.
"""
from decimal import Decimal as D, getcontext
from pathlib import Path
import numpy as np,json,time
getcontext().prec=110
root=Path(__file__).resolve().parents[1];work=root/'work'
ONE=D(1); EPS=D('1e-95');CUT=D('1e-115')
def atan_inv(n):
 x=ONE/D(n); xx=-x*x; power=x;s=power;k=0
 while True:
  k+=1;power*=xx;term=power/D(2*k+1);s+=term
  if abs(term)<CUT: return s
pi=16*atan_inv(5)-4*atan_inv(239)
def sc(x):
 x2=-x*x;ts=x;tc=ONE;s=ts;c=tc;k=0
 while True:
  k+=1;ts*=x2/D((2*k)*(2*k+1));tc*=x2/D((2*k-1)*(2*k));s+=ts;c+=tc
  if abs(ts)+abs(tc)<CUT:return s,c
N=4096;t=time.time();sd=[];cd=[]
for j in range(N+1):
 s,c=sc(2*pi*D(j)/D(N));sd.append(s);cd.append(c)
def bounds(seq):
 lo=np.array([float(x-EPS) for x in seq]);hi=np.array([float(x+EPS) for x in seq]);return np.nextafter(lo,-np.inf),np.nextafter(hi,np.inf)
sl,sh=bounds(sd);cl,ch=bounds(cd)
np.savez(work/'independent_trig4096.npz',cl=cl,ch=ch,sl=sl,sh=sh)
# Compare values as exact Decimal representations of binary64 endpoints.
orig=np.load(root/'data/trig_cc_exact.npz');old4=np.load(work/'trig4096_mp_iv.npz')
results={}
for prefix,n in [('c1024',1024),('s1024',1024),('c2048',2048),('s2048',2048)]:
 vals=cd if prefix[0]=='c' else sd;bad=[]
 for j in range(n):
  val=vals[j*(N//n)];lo=D.from_float(float(orig[prefix+'lo'][j]));hi=D.from_float(float(orig[prefix+'hi'][j]));
  if val+EPS<lo or val-EPS>hi:bad.append(j)
 results[prefix]={'outside_count':len(bad),'first_indices':bad[:12]}
for pref,vals in [('c',cd),('s',sd)]:
 bad=[];maxdist=D(0)
 for j,val in enumerate(vals):
  lo=D.from_float(float(old4[pref+'l'][j]));hi=D.from_float(float(old4[pref+'h'][j]));dist=max(lo-val,val-hi,D(0))
  if dist>EPS:bad.append(j);maxdist=max(maxdist,dist)
 results['generated4096_'+pref]={'outside_count':len(bad),'first_indices':bad[:12],'max_violation':str(maxdist)}
for k,v in results.items():print(k,v,flush=True)
# Exact rational cosine formula, evaluated in high-precision Decimal.
wd=[];end=ONE/D(2048*2048-1);wd.append(end)
for j in range(1,2048):
 sm=sum((2*cd[(2*j*k)%4096]/D(4*k*k-1) for k in range(1,1024)),D(0))
 wd.append((ONE-sm-((-1)**j)*end)/D(1024))
wd.append(end);wl,wh=bounds(wd)
np.savez(work/'independent_cc2048.npz',lo=wl,hi=wh)
json.dump({'pi':str(pi),'z':[str(x) for x in cd[:2049]],'st':[str(x) for x in sd[:2049]],'epsilon':str(EPS)},open(work/'independent_decimal_nodes.json','w'))
json.dump(results,open(work/'independent_table_audit.json','w'),indent=2)
print('weights sum',sum(wd),'seconds',time.time()-t,flush=True)
