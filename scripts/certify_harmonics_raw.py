import json,numpy as np,math,os
from fractions import Fraction as F
HERE=os.path.dirname(os.path.abspath(__file__))
ROOT=os.path.dirname(HERE)
WORK=os.path.join(ROOT,'work')
os.makedirs(WORK,exist_ok=True)
D=json.load(open(os.path.join(ROOT,'data','werner_rational_v6910_candidate_deg111.json'))); weights=[F(s) for s in D['weights']];wf=np.array([float(x) for x in weights]);wr=np.array([abs(float(F.from_float(float(x))-x))+np.spacing(abs(float(x)) if x else 1.0) for x in weights])
T2=np.load(os.path.join(ROOT,'data','trig_cc_exact.npz'));T4=np.load(os.path.join(WORK,'trig4096_mp_iv.npz'));W=np.load(os.path.join(WORK,'cc2048_weights_iv.npz'))
LMAX=151;Nz=2049;R=42;NOR=21
Fmid=np.memmap(os.path.join(WORK,'Fmid_cc2048_direct.dat'),mode='r',dtype=np.complex128,shape=(LMAX+1,Nz,R));META=np.load(os.path.join(WORK,'cc2048_direct_meta.npz'));TOTAL=META['total'];z=META['z'];st=META['st'];wmid=META['wmid'];wlo=W['lo'];whi=W['hi'];wrad=(whi-wlo)/2
zrad=(T4['ch'][:Nz]-T4['cl'][:Nz])/2; sinlo=np.maximum(0,T4['sl'][:Nz]);
pilo=float(T2['pilo'][0]);pihi=float(T2['pihi'][0]);pimid=(pilo+pihi)/2
# previously established bound: degree-9 reconstruction interval <7.49e-12; derivative of response wrt g <630.606.
# Add >2e-9 for all central interpolation/response arithmetic. Conservative uniform response error:
ETA_F=2.0e-8
# The direct preceding stage derives a much smaller response enclosure.
assert float(META['eta_f']) < ETA_F, (float(META['eta_f']), ETA_F)
# exact phi DFT perturbation from response interval, plus a conservative 1e-10 for binary64 FFT and dphi midpoint arithmetic.
EF=2*pihi*ETA_F+1.0e-9
# total energy perturbation: exact-degree CC quadrature, pointwise |delta(f^2)|<=2eta+eta^2 over area 4pi, plus generous numerical/weight padding.
ETOT=4*pihi*(2*ETA_F+ETA_F*ETA_F)+2e-9
# denominator 1/sqrt(2pi) and 1/(4pi)
s2lo=math.sqrt(2*pilo);s2hi=math.sqrt(2*pihi); invs2_mid=1/math.sqrt(2*pimid); invs2_err=max(abs(invs2_mid-1/s2lo),abs(invs2_mid-1/s2hi))+2e-16
fourpi_lo=4*pilo;fourpi_hi=4*pihi
ells=list(range(1,LMAX+1,2));li={l:i for i,l in enumerate(ells)};nL=len(ells)
Elo=np.zeros((nL,R));Oc=np.zeros((nL,NOR));Oe=np.zeros((nL,NOR))
# normalized associated Legendre diagonals at central z
Pdiag=[None]*(LMAX+1);Pdiag[0]=np.full(Nz,1/math.sqrt(2.0))
for m in range(1,LMAX+1):Pdiag[m]=-math.sqrt(1+1/(2*m))*st*Pdiag[m-1]
U=2**-53
for m in range(LMAX+1):
 pprev=None;pcur=Pdiag[m].copy();l=m
 while l<=LMAX:
  if l>=1 and l%2==1:
   # central z projection for all responses
   Fm=np.asarray(Fmid[m,:,:]) # Nz x R
   base=wmid*pcur
   Nmid=np.sum(base[:,None]*Fm,axis=0)
   h=Nmid*invs2_mid
   # Pbar node error: exact-node displacement via the same derivative envelope as original verifier + recurrence arithmetic cushion.
   dU=math.sqrt(l*(l+1)*(2*l+1)/2) if l else 0.0
   ep=np.full(Nz,2e-11*(l+1))
   mask=sinlo>0
   if l: ep[mask]+=dU*zrad[mask]/sinlo[mask]
   # numerator complex-magnitude error
   Fabs=np.abs(Fm); Pabs=np.abs(pcur)
   termF=np.sum(np.abs(wmid)*Pabs)*EF
   termP=np.sum(np.abs(wmid)[:,None]*ep[:,None]*(Fabs+EF),axis=0)
   termW=np.sum(wrad[:,None]*(Pabs[:,None]+ep[:,None])*(Fabs+EF),axis=0)
   # dot roundoff, very conservative 10*gamma_n absolute-sum envelope
   gam=(2*Nz*U)/(1-2*Nz*U)
   doterr=10*gam*np.sum(np.abs(base)[:,None]*Fabs,axis=0)+2e-12
   Nerr=termF+termP+termW+doterr
   herr=(Nerr*(1/s2lo))+invs2_err*(np.abs(Nmid)+Nerr)+2e-12
   idx=li[l]; fac=1 if m==0 else 2
   mag=np.abs(h); low=np.maximum(0,mag-herr);Elo[idx]+=fac*low*low
   a=h[:NOR];b=h[NOR:];da=herr[:NOR];db=herr[NOR:]
   Oc[idx]+=fac*np.real(a*np.conj(b));Oe[idx]+=fac*(np.abs(a)*db+np.abs(b)*da+da*db)
  # advance recurrence
  if l==LMAX: break
  if l==m:
   pnext=z*math.sqrt(2*m+3)*pcur
  else:
   ll=l+1
   A=math.sqrt((4*ll*ll-1)/(ll*ll-m*m));B=math.sqrt(((2*ll+1)*((ll-1)*(ll-1)-m*m))/((2*ll-3)*(ll*ll-m*m)))
   pnext=A*z*pcur-B*pprev
  pprev,pcur=pcur,pnext;l+=1
 if m%20==0: print('m',m,flush=True)
# c_l mixture intervals from orbit overlap intervals /4pi
clo=[];chi=[];cc=[];ce=[]
for idx,l in enumerate(ells):
 oc=Oc[idx];oe=Oe[idx]
 # per-orbit correlation central and error including 4pi interval
 c0=oc/(4*pimid)
 # division denominator perturb + overlap interval
 denomerr=max(abs(1/fourpi_lo-1/(4*pimid)),abs(1/fourpi_hi-1/(4*pimid)))
 er=oe/fourpi_lo + np.abs(oc)*denomerr + 2e-15
 mixc=float(np.dot(wf,c0));mixe=float(np.dot(np.abs(wf),er)+np.dot(wr,np.abs(c0)+er)+2e-15)
 clo.append(mixc-mixe);chi.append(mixc+mixe);cc.append(mixc);ce.append(mixe)
 if l%10==1: print('l',l,'c',mixc,'err',mixe,'EloMinMax',Elo[idx].min(),Elo[idx].max(),flush=True)
clo=np.array(clo);chi=np.array(chi);cc=np.array(cc);ce=np.array(ce)
# resolved residual upper
Dres=0.0
for idx,l in enumerate(ells):
 if l>=3:Dres+=max(abs(clo[idx]),abs(chi[idx]))
# tail from total upper - cumulative low energy lower
cumlo=np.sum(Elo,axis=0);tot_hi=TOTAL+ETOT;tailE=np.maximum(0,tot_hi-cumlo)
tail=0.0
for r in range(NOR):
 term=math.sqrt(tailE[r]*tailE[NOR+r])/fourpi_lo
 # exact positive weight interval upper
 wup=wf[r]+wr[r];tail+=wup*term
# v0 interval and corrected lower
v0lo=-chi[0];v0hi=-clo[0];Dtot=Dres+tail;K=1.783;vc=v0lo/(1+K*Dtot)
out={'eta_f':ETA_F,'EF':EF,'ETOT_rad':ETOT,'c1_center':cc[0],'c1_error':ce[0],'v0_lower':v0lo,'v0_upper':v0hi,'resolved_D_upper':Dres,'tail_upper':tail,'D_total_upper':Dtot,'v_cert_lower':vc,'LMAX':LMAX,'tailE':tailE.tolist(),'coeff_intervals':[{'ell':l,'lo':float(clo[i]),'hi':float(chi[i]),'center':float(cc[i]),'radius':float(ce[i])} for i,l in enumerate(ells)]}
json.dump(out,open(os.path.join(WORK,'werner_harmonic_certificate_deg111_regenerated.json'),'w'),indent=2)
print('TOTAL central minmax',TOTAL.min(),TOTAL.max(),'ETOT',ETOT)
print('c1',clo[0],chi[0],'v0',v0lo,v0hi)
print('Dres',Dres,'tail',tail,'Dtot',Dtot,'VCERT_LOWER',vc)
print('max tailE',tailE.max(),'min',tailE.min())
if vc>0.6910: print('*** CERTIFICATE CLEARS 0.6910 ***')
else: print('*** DOES NOT CLEAR 0.6910 ***')
