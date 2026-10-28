import json,numpy as np,math,time,os
from fractions import Fraction as F
HERE=os.path.dirname(os.path.abspath(__file__))
ROOT=os.path.dirname(HERE)
WORK=os.path.join(ROOT,'work')
os.makedirs(WORK,exist_ok=True)
def get_num_threads(): return 1
D=json.load(open(os.path.join(ROOT,'data','werner_rational_v6910_candidate_deg111.json')));E=[tuple(e) for e in D['exponents']]
rows=D['A']+D['B'];C=np.array([[float(F(s)) for s in row] for row in rows],float);Cr=np.zeros_like(C)
for i,row in enumerate(rows):
 for j,s in enumerate(row):
  q=F(s);x=float(q);Cr[i,j]=abs(float(F.from_float(x)-q))+np.spacing(abs(x) if x else 1.0)
pcQ=[F(s) for s in D['response']['cheb_coeffs']];pc=np.array([float(x) for x in pcQ]);odd=pc[1::2].copy();Lip_u=float(sum(abs(a)*k*k for k,a in enumerate(pcQ)));Lip_g=Lip_u/0.9
T4=np.load(os.path.join(WORK,'trig4096_mp_iv.npz'));T2=np.load(os.path.join(ROOT,'data','trig_cc_exact.npz'));W=np.load(os.path.join(WORK,'cc2048_weights_iv.npz'))
Nz=2049;R=42;Nphi=2048;zlo,zhi=T4['cl'][:Nz],T4['ch'][:Nz];rlo,rhi=T4['sl'][:Nz],T4['sh'][:Nz];z=(zlo+zhi)/2;st=(rlo+rhi)/2
cL,cH,sL,sH=T2['c2048lo'],T2['c2048hi'],T2['s2048lo'],T2['s2048hi'];cm=(cL+cH)/2;sm=(sL+sH)/2;cr=(cH-cL)/2;sr=(sH-sL)/2
idx32=np.arange(0,Nphi,64); U=2**-53
# DFT32 matrices m=0..9, samples k=0..31
kk=np.arange(32);E32=np.empty((10,32),complex);ER32=np.empty((10,32))
for m in range(10):
 ix=(m*idx32)%Nphi;E32[m]=cm[ix]-1j*sm[ix];ER32[m]=np.hypot(cr[ix],sr[ix])
# full reconstruction positive exponential m=1..9
jj=np.arange(Nphi);Efull=np.empty((9,Nphi),complex);ERfull=np.empty((9,Nphi))
for m in range(1,10):
 ix=(m*jj)%Nphi;Efull[m-1]=cm[ix]+1j*sm[ix];ERfull[m-1]=np.hypot(cr[ix],sr[ix])
max_tw_full=float(ERfull.max())
wmid=(W['lo']+W['hi'])/2;pimid=(float(T2['pilo'][0])+float(T2['pihi'][0]))/2;dphi=2*pimid/Nphi
LMAX=151;Fmid=np.memmap(os.path.join(WORK,'Fmid_cc2048_direct.dat'),mode='w+',dtype=np.complex128,shape=(LMAX+1,Nz,R));Ez=np.memmap(os.path.join(WORK,'Ez_cc2048_direct.dat'),mode='w+',dtype=np.float64,shape=(Nz,R))
def evalodd(g,co):
 u=g/0.9
 x=2*u*u-1.0
 r0=u
 s=co[0]*r0
 r1=u*(2*x-1.0)
 s+=co[1]*r1
 for k in range(2,len(co)):
  r2=2*x*r1-r0
  s+=co[k]*r2
  r0,r1=r1,r2
 return s
_=evalodd(np.zeros((2,3)),odd);print('threads',get_num_threads(),'Lip_g',Lip_g,flush=True)
INF=np.inf
def down(x):return np.nextafter(x,-INF)
def up(x):return np.nextafter(x,INF)
def imul(a,b):
 p0=a[0]*b[0];p1=a[0]*b[1];p2=a[1]*b[0];p3=a[1]*b[1];return down(np.minimum(np.minimum(p0,p1),np.minimum(p2,p3))),up(np.maximum(np.maximum(p0,p1),np.maximum(p2,p3)))
def imidrad(a):
 m=a[0]/2+a[1]/2;r=up(np.maximum(m-a[0],a[1]-m));return m,r
max_samp_rad=0.;max_coef_rad=0.;max_g_rad=0.;t0=time.time();bs=16
for j0 in range(0,Nz,bs):
 inds=np.arange(j0,min(Nz,j0+bs));bz=len(inds);# interval small-grid coords
 CP=(cL[idx32][None,:],cH[idx32][None,:]);SP=(sL[idx32][None,:],sH[idx32][None,:]);RR=(rlo[inds,None],rhi[inds,None]);ZZ=(zlo[inds,None],zhi[inds,None]);X=imul(RR,CP);Y=imul(RR,SP);Zb=(np.broadcast_to(ZZ[0],(bz,32)),np.broadcast_to(ZZ[1],(bz,32)))
 def powers(A):
  o=[(np.ones_like(A[0]),np.ones_like(A[0])),A]
  for _ in range(2,10):o.append(imul(o[-1],A))
  return o
 XP,YP,ZP=powers(X),powers(Y),powers(Zb);P=bz*32;Mm=np.empty((P,55));Mr=np.empty_like(Mm)
 for k,(a,b,c) in enumerate(E):
  I=imul(imul(XP[a],YP[b]),ZP[c]);mm,rr=imidrad(I);Mm[:,k]=mm.ravel();Mr[:,k]=rr.ravel()
 gm=(Mm@C.T).reshape(bz,32,R);rad=(Mr@np.abs(C).T+(np.abs(Mm)+Mr)@Cr.T).reshape(bz,32,R);gam=(2*55*U)/(1-2*55*U);rad+= (2*gam*(np.abs(Mm)@np.abs(C).T+1e-300)).reshape(bz,32,R);max_samp_rad=max(max_samp_rad,float(rad.max()))
 # exact DFT coefficients intervals around direct central coefficient computation
 Ac=np.einsum('mk,bkr->bmr',E32,gm,optimize=True)/32.0
 # interval/error for each complex coefficient
 Ar=np.empty((bz,10,R))
 for m in range(10):
  # input interval + twiddle intervals + product cross + dot roundoff
  twr=ER32[m][None,:,None];twabs=np.abs(E32[m])[None,:,None]
  er=np.sum(twabs*rad+np.abs(gm)*twr+rad*twr,axis=1)/32.0
  gamma=(2*32*U)/(1-2*32*U);er+=10*gamma*np.sum(np.abs(gm)*twabs,axis=1)/32.0+2e-14
  Ar[:,m,:]=er
 max_coef_rad=max(max_coef_rad,float(Ar.max()))
 # direct central reconstruction on full phi
 G=Ac[:,0,:].real[:,None,:]+2*np.real(np.einsum('bmr,mk->bkr',Ac[:,1:,:],Efull,optimize=True))
 # rigorous uniform reconstruction radius per z,response (independent of phi using max twiddle rad)
 gr=Ar[:,0,:].copy()
 gr+=2*np.sum(Ar[:,1:,:]+np.abs(Ac[:,1:,:])*max_tw_full+Ar[:,1:,:]*max_tw_full,axis=1)
 gamma=(2*9*U)/(1-2*9*U);gr+=20*gamma*np.sum(np.abs(Ac[:,1:,:]),axis=1)+2e-13
 max_g_rad=max(max_g_rad,float(gr.max()))
 Q=evalodd(G,odd);Ez[inds,:]=dphi*np.sum(Q*Q,axis=1);fq=np.fft.fft(Q,axis=1)*dphi;Fmid[:,inds,:]=np.transpose(fq[:,:LMAX+1,:],(1,0,2))
 if j0%256==0:print('z',j0,'elapsed',time.time()-t0,'srad',max_samp_rad,'Ar',max_coef_rad,'Grad',max_g_rad,'Gmax',np.max(abs(G)),flush=True)
Fmid.flush();Ez.flush();TOTAL=np.sum(wmid[:,None]*np.asarray(Ez),axis=0)
# response arithmetic bound inflated beyond degree*roundoff; exact coefficient rounding included separately below.
pc_round=sum(abs(float(F.from_float(float(x))-q)) for x,q in zip(pc,pcQ))
RESP_ARITH=2e-11
ETA_F=Lip_g*max_g_rad+pc_round+RESP_ARITH
np.savez(os.path.join(WORK,'cc2048_direct_meta.npz'),total=TOTAL,z=z,st=st,wmid=wmid,dphi=dphi,max_sample_rad=max_samp_rad,max_fourier_coeff_rad=max_coef_rad,max_g_rad=max_g_rad,Lip_g=Lip_g,pc_round=pc_round,response_arith=RESP_ARITH,eta_f=ETA_F)
print('DONE',time.time()-t0,'maxsamp',max_samp_rad,'maxcoef',max_coef_rad,'maxgrad',max_g_rad,'pcround',pc_round,'ETA_F',ETA_F,'TOTAL',TOTAL.min(),TOTAL.max(),flush=True)
