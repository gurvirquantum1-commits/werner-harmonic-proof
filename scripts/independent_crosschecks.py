from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,getcontext
import json,numpy as np,math
root=Path(__file__).resolve().parents[1];w=root/'work';data=json.load(open(root/'data/werner_rational_v6910_candidate_deg111.json'));getcontext().prec=100
coq=[F(s) for s in data['response']['cheb_coeffs']];cod=[D(x.numerator)/D(x.denominator) for x in coq];co=np.array([float(x) for x in coq]);u=np.linspace(-.973,.973,3001);xx=2*u*u-1;r0=u;q=co[1]*r0;r1=u*(2*xx-1);q+=co[3]*r1
for k,c in enumerate(co[5::2],2):r2=2*xx*r1-r0;q+=c*r2;r0,r1=r1,r2
md=0.
for uf,qf in zip(u,q):
 ud=D.from_float(float(uf));b1=D(0);b2=D(0)
 for a in cod[:0:-1]:b0=2*ud*b1-b2+a;b2,b1=b1,b0
 exact=ud*b1-b2+cod[0];md=max(md,float(abs(exact-D.from_float(float(qf)))))
T=np.load(w/'independent_trig4096.npz');W=np.load(w/'independent_cc2048.npz');wm=(W['lo']+W['hi'])/2;z=(T['cl'][:2049]+T['ch'][:2049])/2;st=(T['sl'][:2049]+T['sh'][:2049])/2;fm=np.memmap(w/'independent_Fmid.dat',mode='r',dtype=np.complex128,shape=(152,2049,42));moment=np.array([(wm*st)@fm[1].real,-(wm*st)@fm[1].imag,(wm*z)@fm[0].real]);weights=np.array([float(F(x)) for x in data['weights']]);cart=3/(16*math.pi**2)*np.sum(weights*np.sum(moment[:,:21]*moment[:,21:],axis=0));cert=json.load(open(w/'independent_certificate.json'))
r={'3001_response_samples_max_100digit_error':md,'dipole_via_cartesian_moments':cart,'dipole_via_harmonics':cert['c1_center'],'dipole_discrepancy':abs(cart-cert['c1_center'])}
print(json.dumps(r,indent=2));assert md<2e-11 and abs(cart-cert['c1_center'])<1e-13
json.dump(r,open(w/'independent_crosschecks.json','w'),indent=2)
