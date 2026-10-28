"""Projection audit using independent direct DFT and Decimal Legendre balls.
Quadrature is exact: phi frequencies for response*Y <=1150, for response^2
<=1998, both <2048. After azimuthal integration polar polynomials have
respective degrees <=1150 and <=1998, both <=CC degree2048.
Additional 1e-9 total-energy/resolved-energy guards and 1e-11 orbit overlap
guards dominate gamma_12000 roundoff of positive energy and overlap sums
(bounded by ~13), so their subtraction and mixture need no exact cancellation.
"""
from pathlib import Path
from fractions import Fraction as Q
import numpy as np,math,json
root=Path(__file__).resolve().parents[1];work=root/'work';D=json.load(open(root/'data/werner_rational_v6910_candidate_deg111.json'));nodes=json.load(open(work/'independent_decimal_nodes.json'));pi=float(nodes['pi']);plo=np.nextafter(pi,-np.inf);phi=np.nextafter(pi,np.inf);eta=2e-8;EF=2*phi*eta+1e-9;ETOT=4*phi*(2*eta+eta*eta)+1e-9
W=np.load(work/'independent_cc2048.npz');wm=(W['lo']+W['hi'])/2;wr=(W['hi']-W['lo'])/2+2**-53*abs(wm)
P=np.memmap(work/'independent_Pbar.dat',mode='r',dtype=np.float64,shape=(152,152,2049));Pr=np.load(work/'independent_Pbar_radius.npy');Fm=np.memmap(work/'independent_Fmid.dat',mode='r',dtype=np.complex128,shape=(152,2049,42));meta=np.load(work/'independent_meta.npz');weights=np.array([float(Q(s)) for s in D['weights']]);weights_err=2**-52*weights
ells=list(range(1,152,2));Elo=np.zeros((76,42));Oc=np.zeros((76,21));Oe=np.zeros((76,21));inv=1/math.sqrt(2*pi)
# Check normalized Legendre recurrence from supplied regeneration against
# independent Decimal true-node projections, including its claimed envelope.
oldT=np.load(work/'trig4096_mp_iv.npz');z=(oldT['cl'][:2049]+oldT['ch'][:2049])/2;st=(oldT['sl'][:2049]+oldT['sh'][:2049])/2;old_diag=np.full(2049,1/math.sqrt(2));max_ratio=0.;maxdiff=0.
for m in range(152):
 if m:old_diag=-math.sqrt(1+1/(2*m))*st*old_diag
 oldp=old_diag.copy();oldprev=None
 for l in range(m,152):
  p=P[m,l];pr=float(Pr[m,l]);delta=abs(oldp-p)+pr;old_allowance=2e-11*(l+1);max_ratio=max(max_ratio,float(delta.max()/old_allowance));maxdiff=max(maxdiff,float(delta.max()))
  if l%2:
   j=(l-1)//2;fac=1 if m==0 else 2;base=wm*p;F=np.asarray(Fm[m]);N=np.sum(base[:,None]*F,axis=0);h=N*inv
   # Bounds for all projection ingredients plus 2049-term complex dot.
   herr=(np.sum(abs(wm)*abs(p))*EF + np.sum(abs(wm))*pr*(np.max(abs(F),axis=0)+EF) + np.sum(wr*(abs(p)+pr))*(np.max(abs(F),axis=0)+EF)+1e-10)*inv+1e-12
   assert np.max(abs(h))<4
   Elo[j]+=fac*np.maximum(0,abs(h)-herr)**2
   a,b=h[:21],h[21:];da,db=herr[:21],herr[21:]
   Oc[j]+=fac*np.real(a*np.conj(b));Oe[j]+=fac*(abs(a)*db+abs(b)*da+da*db)
  if l==151:break
  if l==m:nxt=z*math.sqrt(2*m+3)*oldp
  else:
   n=l+1;A=math.sqrt((4*n*n-1)/(n*n-m*m));B=math.sqrt(((2*n+1)*((n-1)**2-m*m))/((2*n-3)*(n*n-m*m)));nxt=A*z*oldp-B*oldprev
  oldprev,oldp=oldp,nxt
 if m%30==0:print('m',m,flush=True)
cc=[];ce=[]
for j,l in enumerate(ells):
 orbit=Oc[j]/(4*pi);err=(Oe[j]+1e-11)/(4*plo)+1e-14
 cc.append(float(weights@orbit));ce.append(float(weights@err+weights_err@(abs(orbit)+err)+1e-13))
cc=np.array(cc);ce=np.array(ce);lo=cc-ce;hi=cc+ce
retained=np.sum(Elo,axis=0)-1e-9;total_hi=meta['total']+ETOT;tailE=np.maximum(0,total_hi-retained)
tail=float(np.sum((weights+weights_err)*np.sqrt(tailE[:21]*tailE[21:]))/(4*plo))+1e-12
resolved=float(np.sum(np.maximum(abs(lo[1:]),abs(hi[1:]))))+1e-12
v0=-hi[0];dtot=resolved+tail;v=v0/(1+1.783*dtot)
original=json.load(open(root/'data/werner_harmonic_certificate_deg111_final.json'));center_discrepancy=max(abs(cc[j]-x['center']) for j,x in enumerate(original['coeff_intervals']));inside=all(x['lo']<=cc[j]<=x['hi'] for j,x in enumerate(original['coeff_intervals']))
contained=all(o['lo']<=lo[j]<=hi[j]<=o['hi'] for j,o in enumerate(original['coeff_intervals']))
tails_contained=bool(np.all(tailE<=np.array(original['tailE'])))
assert contained and tails_contained
result={'all_independent_intervals_inside_archived_intervals':contained,'all_independent_tail_bounds_below_archived_bounds':tails_contained,'minimum_archived_tail_energy_slack':float(np.min(np.array(original['tailE'])-tailE)),'method':'independent Decimal trig/CC/Legendre; direct polynomial evaluation; direct 2048-term DFT','v0_lower':v0,'c1_center':cc[0],'c1_radius':ce[0],'resolved_D_upper':resolved,'tail_upper':tail,'D_total_upper':dtot,'v_cert_lower':v,'max_supplied_Legendre_error':maxdiff,'max_fraction_of_supplied_Legendre_allowance':max_ratio,'max_original_center_discrepancy':center_discrepancy,'independent_centers_in_original_intervals':inside,'tailE':tailE.tolist(),'coeff_intervals':[{'ell':l,'center':float(cc[j]),'radius':float(ce[j]),'lo':float(lo[j]),'hi':float(hi[j])} for j,l in enumerate(ells)],'error_constants':{'response':eta,'direct_DFT':EF,'energy':ETOT,'energy_sum_guard':1e-9,'orbit_overlap_guard':1e-11}}
json.dump(result,open(work/'independent_certificate.json','w'),indent=2)
for k,vv in result.items():
 if k not in ['tailE','coeff_intervals']:print(k,vv)
assert v0>0.6914213 and dtot<0.0001213
bound=Q(6914213,10_000_000)/(1+Q(1783,1000)*Q(1213,10_000_000));assert bound==Q(6914213000,10002162779)
print('PASS: independently reconstructed conservative bounds imply exact manuscript bound',bound,float(bound))
