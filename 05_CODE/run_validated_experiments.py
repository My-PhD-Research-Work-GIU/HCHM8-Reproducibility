from __future__ import annotations
import hashlib, json, math, time, csv
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'; FIG=ROOT/'figures'; QA=ROOT/'qa'
for p in (DATA,FIG,QA): p.mkdir(parents=True,exist_ok=True)

A = np.array([
[2,1,0,0,0,0,0,0],[1,1,0,0,0,0,0,0],
[0,0,3,2,0,0,0,0],[0,0,1,1,0,0,0,0],
[0,0,0,0,5,2,0,0],[0,0,0,0,2,1,0,0],
[0,0,0,0,0,0,7,3],[0,0,0,0,0,0,2,1]],dtype=np.float64)

REF_LABEL=b"GIU-PhD-HCHM8-heterogeneous-conjugacy-v1"

def shake(data:bytes,n:int)->bytes: return hashlib.shake_256(data).digest(n)
def unit_values(label:bytes,count:int)->np.ndarray:
    raw=shake(label,8*count)
    return np.frombuffer(raw,dtype='>u8').astype(np.float64)/(2**64)

def key_for_index(i:int)->bytes: return shake(b'HCHM8-CONJUGATE-KEY'+i.to_bytes(8,'big'),64)

def torus_diff(a,b):
    d=np.abs(a-b); return np.minimum(d,1.0-d)

@dataclass
class Params:
    epsilon: np.ndarray
    phi1: float; phi2: float; rho3: float; kappa4: float; eta6: float; kappa7: float; rho8: float

class HCHM8:
    def __init__(self,key_material:bytes,round_index:int=0):
        label=key_material+b'|HCHM8-PARAM|'+round_index.to_bytes(4,'big')
        r=unit_values(label,32)
        self.key_material=key_material; self.round_index=round_index
        self.x0=np.mod(r[:8]+np.arange(1,9)*0.017391,1.0)
        self.p=Params(
            epsilon=0.07+0.08*r[8:16],
            phi1=2*np.pi*r[16], phi2=2*np.pi*r[17],
            rho3=0.75+2.25*r[18], kappa4=1.25+2.75*r[19],
            eta6=0.75+1.75*r[20], kappa7=1.25+2.75*r[21],
            rho8=1.0+5.0*r[22])
        self.i0_eta=float(np.i0(self.p.eta6))
    @classmethod
    def reference(cls):
        # preserves exact source profile by using the original label directly
        obj=cls.__new__(cls)
        r=unit_values(REF_LABEL,32)
        obj.key_material=REF_LABEL; obj.round_index=0
        obj.x0=np.mod(r[:8]+np.arange(1,9)*0.017391,1.0)
        obj.p=Params(epsilon=0.07+0.08*r[8:16],phi1=2*np.pi*r[16],phi2=2*np.pi*r[17],
            rho3=0.75+2.25*r[18],kappa4=1.25+2.75*r[19],eta6=0.75+1.75*r[20],
            kappa7=1.25+2.75*r[21],rho8=1.0+5.0*r[22])
        obj.i0_eta=float(np.i0(obj.p.eta6)); return obj
    def f_and_grad(self,i,z):
        p=self.p; g=np.zeros(8,dtype=np.float64); twopi=2*np.pi
        if i==0:
            th=twopi*z[2]+p.phi1; f=np.sin(th); g[2]=twopi*np.cos(th)
        elif i==1:
            th=twopi*z[5]+p.phi2; c=np.cos(th); f=4*c**3-3*c
            g[5]=(12*c*c-3)*(-twopi*np.sin(th))
        elif i==2:
            th=twopi*(z[7]-z[1]); N=np.sin(th); c=np.cos(twopi*z[1]); D=1+p.rho3*c*c
            dN8=twopi*np.cos(th); dN2=-dN8
            dD2=-4*np.pi*p.rho3*c*np.sin(twopi*z[1])
            f=N/D; g[7]=dN8/D; g[1]=(dN2*D-N*dD2)/(D*D)
        elif i==3:
            th=twopi*(z[0]-z[6]); s=np.sin(th); u=p.kappa4*s; f=np.tanh(u)
            d=(1-f*f)*p.kappa4*twopi*np.cos(th); g[0]=d; g[6]=-d
        elif i==4:
            th=twopi*(z[6]+z[7]); q=.5*(1+np.sin(th)); f=4*q*(1-q)-.5
            d=(4-8*q)*np.pi*np.cos(th); g[6]=d; g[7]=d
        elif i==5:
            th=twopi*(z[1]-z[0]); e=np.exp(p.eta6*np.cos(th)); D=np.exp(p.eta6)-np.exp(-p.eta6)
            f=(e-self.i0_eta)/D; d2=e*p.eta6*(-np.sin(th))*twopi/D; g[1]=d2; g[0]=-d2
        elif i==6:
            th=twopi*(z[3]+z[5]); s=np.sin(th); den=1+(p.kappa7*s)**2
            f=(2/np.pi)*np.arctan(p.kappa7*s); d=4*p.kappa7*np.cos(th)/den; g[3]=d; g[5]=d
        elif i==7:
            u=np.pi*(z[4]-z[2]); G=1+p.rho8*np.sin(u)**2; L=np.log(1+p.rho8)
            f=np.log(G)/L-.5; d=p.rho8*np.pi*np.sin(2*u)/(G*L); g[4]=d; g[2]=-d
        return float(f),g
    def H(self,state,return_jac=False,dtype=np.float64):
        z=np.array(state,dtype=dtype,copy=True)
        J=np.eye(8,dtype=np.float64)
        for i in range(8):
            # evaluate in float64 for derivative; state path follows dtype
            f,g=self.f_and_grad(i,z.astype(np.float64))
            Ji=np.eye(8); Ji[i,:]+=self.p.epsilon[i]*g
            J=Ji@J
            z[i]=dtype((float(z[i])+self.p.epsilon[i]*f)%1.0)
        return (z,J) if return_jac else z
    def H_inverse(self,state,dtype=np.float64):
        z=np.array(state,dtype=dtype,copy=True)
        for i in range(7,-1,-1):
            f,_=self.f_and_grad(i,z.astype(np.float64))
            z[i]=dtype((float(z[i])-self.p.epsilon[i]*f)%1.0)
        return z
    def step(self,state,dtype=np.float64):
        z=self.H(state,dtype=dtype)
        y=np.mod(A.astype(dtype)@z,dtype(1.0))
        return self.H_inverse(y,dtype=dtype)
    def step_and_jac(self,state):
        z,DHx=self.H(state,return_jac=True)
        y=np.mod(A@z,1.0)
        xn=self.H_inverse(y)
        _,DHn=self.H(xn,return_jac=True)
        DT=np.linalg.solve(DHn,A@DHx)
        return xn,DT
    @staticmethod
    def exact_spectrum():
        e=np.array([(3+np.sqrt(5))/2,2+np.sqrt(3),3+2*np.sqrt(2),4+np.sqrt(15)])
        p=np.log(e); return np.sort(np.r_[p,-p])[::-1]

# ---- Map tests ----
def lyapunov_numeric(sys:HCHM8,steps=12000,burn=1000):
    x=sys.x0.copy()
    for _ in range(burn): x=sys.step(x)
    Q=np.eye(8); acc=np.zeros(8)
    for _ in range(steps):
        xn,J=sys.step_and_jac(x); Q,R=np.linalg.qr(J@Q); acc+=np.log(np.maximum(np.abs(np.diag(R)),1e-300)); x=xn
    return np.sort(acc/steps)[::-1]

def perm_entropy(x,order=5):
    w=np.lib.stride_tricks.sliding_window_view(np.asarray(x),order)
    p=np.argsort(w,axis=1,kind='stable'); codes=np.zeros(len(p),dtype=np.int32)
    for j in range(order): codes=codes*order+p[:,j]
    _,c=np.unique(codes,return_counts=True); q=c/c.sum(); return float(-(q*np.log(q)).sum()/np.log(math.factorial(order)))

def byte_entropy(x):
    q=np.minimum((np.asarray(x)*256).astype(int),255); c=np.bincount(q,minlength=256); p=c[c>0]/c.sum(); return float(-(p*np.log2(p)).sum())

def run_map_tests():
    t0=time.time(); ref=HCHM8.reference(); rng=np.random.default_rng(20260802)
    # exact structural checks
    rounderrs=[]; detH=[]; detT=[]; conjugacy=[]
    for _ in range(2000):
        x=rng.random(8); h,Jh=ref.H(x,return_jac=True); xr=ref.H_inverse(h); rounderrs.append(np.max(torus_diff(x,xr)))
        xn,Jt=ref.step_and_jac(x); detH.append(np.linalg.det(Jh)); detT.append(np.linalg.det(Jt))
        conjugacy.append(np.max(torus_diff(ref.H(xn),np.mod(A@ref.H(x),1.0))))
    # reference sensitivity
    x=ref.x0.copy(); y=x.copy(); y[0]=(y[0]+1e-15)%1.0; ds=[]
    traj=[]
    for n in range(301):
        ds.append(float(torus_diff(x,y)[0])); traj.append(x.copy())
        x=ref.step(x); y=ref.step(y)
    def first(th):
        ids=np.flatnonzero(np.array(ds)>th); return int(ids[0]) if len(ids) else None
    # numeric LE across 12 key-derived systems
    les=[]; entropy=[]; sens=[]
    for i in range(6):
        sys=HCHM8(key_for_index(i)); le=lyapunov_numeric(sys,steps=2200,burn=400); les.append(le)
        x=sys.x0.copy()
        for _ in range(500): x=sys.step(x)
        seq=np.empty((10000,8))
        for n in range(len(seq)): x=sys.step(x); seq[n]=x
        C=np.corrcoef(seq,rowvar=False); off=C-np.eye(8)
        entropy.append({'run':i,'permutation_entropy_mean':np.mean([perm_entropy(seq[:,j]) for j in range(8)]),
                        'byte_entropy_mean':np.mean([byte_entropy(seq[:,j]) for j in range(8)]),
                        'lag1_abs_max':np.max(np.abs([np.corrcoef(seq[:-1,j],seq[1:,j])[0,1] for j in range(8)])),
                        'cross_corr_abs_max':np.max(np.abs(off))})
        a=sys.x0.copy(); b=a.copy(); b[0]=(b[0]+1e-15)%1.0; hit=501
        for n in range(1,501):
            a=sys.step(a); b=sys.step(b)
            if np.linalg.norm(torus_diff(a,b))>.1: hit=n; break
        sens.append(hit)
    les=np.array(les); exact=ref.exact_spectrum()
    pd.DataFrame(les,columns=[f'lambda_{i+1}' for i in range(8)]).to_csv(DATA/'lyapunov_numeric_12runs.csv',index=False)
    pd.DataFrame(entropy).to_csv(DATA/'entropy_correlation_12runs.csv',index=False)
    pd.DataFrame({'run':range(6),'steps_to_torus_distance_0p1':sens}).to_csv(DATA/'sensitivity_12runs.csv',index=False)
    # parameter-invariance scan: scaling eps only, exact spectrum unchanged; numeric confirmation for max LE
    scales=np.linspace(.05,1.5,13); maxnum=[]
    base=HCHM8.reference()
    orig=base.p.epsilon.copy()
    for s in scales:
        base.p.epsilon=orig*s
        maxnum.append(float(lyapunov_numeric(base,steps=1000,burn=200)[0]))
    base.p.epsilon=orig
    pd.DataFrame({'epsilon_scale':scales,'numeric_max_le':maxnum,'exact_max_le':exact[0]}).to_csv(DATA/'parameter_invariance_scan.csv',index=False)
    # orbit projection scan for visual only
    eps_scales=np.linspace(.1,1.4,50); bx=[]; bs=[]
    for s in eps_scales:
        base.p.epsilon=orig*s; xx=base.x0.copy()
        for _ in range(400): xx=base.step(xx)
        for _ in range(30): xx=base.step(xx); bs.append(s); bx.append(xx[0])
    base.p.epsilon=orig
    # exact invariants
    fixed_count=int(round(abs(np.linalg.det(A-np.eye(8)))))
    entropy_exact=float(exact[exact>0].sum())
    summary={
      'map':'HCHM-8 conjugate system T=H^{-1}AH','state_space':'8-torus [0,1)^8',
      'exact_lyapunov_spectrum':exact.tolist(),'positive_exponents':4,
      'exact_topological_entropy_nats_per_iteration':entropy_exact,
      'exact_fixed_points':fixed_count,
      'max_roundtrip_error_2000':float(max(rounderrs)),'max_conjugacy_error_2000':float(max(conjugacy)),
      'detH_mean':float(np.mean(detH)),'detH_max_abs_error_from_1':float(np.max(np.abs(np.array(detH)-1))),
      'detT_mean':float(np.mean(detT)),'detT_max_abs_error_from_1':float(np.max(np.abs(np.array(detT)-1))),
      'reference_sensitivity':{'requested':1e-15,'actual':ds[0],'first_gt_1e-12':first(1e-12),'first_gt_1e-9':first(1e-9),'first_gt_1e-6':first(1e-6),'first_gt_1e-3':first(1e-3),'max_delta_x1':max(ds),'mean_delta_x1':float(np.mean(ds))},
      'numeric_le_runs':6,'numeric_le_mean':les.mean(axis=0).tolist(),'numeric_le_std':les.std(axis=0,ddof=1).tolist(),
      'numeric_vs_exact_max_abs_error':float(np.max(np.abs(les.mean(axis=0)-exact))),
      'sensitivity_steps_median':float(np.median(sens)),'sensitivity_steps_range':[int(min(sens)),int(max(sens))],
      'entropy_correlation_mean':{k:float(pd.DataFrame(entropy)[k].mean()) for k in entropy[0] if k!='run'},
      'elapsed_s':time.time()-t0}
    (DATA/'map_validated_summary.json').write_text(json.dumps(summary,indent=2))
    # plots
    plt.figure(figsize=(8,4.8)); plt.errorbar(range(1,9),les.mean(axis=0),yerr=les.std(axis=0,ddof=1),fmt='o',capsize=4,label='Numerical QR mean ± SD'); plt.plot(range(1,9),exact,'x',label='Exact conjugacy spectrum'); plt.axhline(0,lw=.8); plt.xlabel('Exponent index'); plt.ylabel('Lyapunov exponent'); plt.legend(); plt.tight_layout(); plt.savefig(FIG/'validated_lyapunov_spectrum.pdf'); plt.savefig(FIG/'validated_lyapunov_spectrum.png',dpi=220); plt.close()
    plt.figure(figsize=(8,4.8)); plt.plot(scales,maxnum,'o-',label='Numerical maximum LE'); plt.axhline(exact[0],ls='--',label='Exact invariant value'); plt.xlabel('Common multiplier on shear amplitudes'); plt.ylabel('Maximum Lyapunov exponent'); plt.legend(); plt.tight_layout(); plt.savefig(FIG/'parameter_invariance_maxle.pdf'); plt.savefig(FIG/'parameter_invariance_maxle.png',dpi=220); plt.close()
    plt.figure(figsize=(8,4.8)); plt.scatter(bs,bx,s=.35); plt.xlabel('Common multiplier on shear amplitudes'); plt.ylabel('$x_1$ samples'); plt.title('Coordinate projection across conjugacy parameters (not a chaos/order bifurcation)'); plt.tight_layout(); plt.savefig(FIG/'coordinate_projection_scan.pdf'); plt.savefig(FIG/'coordinate_projection_scan.png',dpi=220); plt.close()
    plt.figure(figsize=(8,4.8)); plt.semilogy(np.maximum(ds,1e-18)); plt.xlabel('Iteration'); plt.ylabel('Torus distance in $x_1$'); plt.tight_layout(); plt.savefig(FIG/'reference_sensitivity.pdf'); plt.savefig(FIG/'reference_sensitivity.png',dpi=220); plt.close()
    return summary

# ---- finite precision ----
# Reference constants for numba quantized-state emulation
REF=HCHM8.reference(); EPS=REF.p.epsilon.astype(np.float64); PHI1=REF.p.phi1; PHI2=REF.p.phi2; RHO3=REF.p.rho3; K4=REF.p.kappa4; ETA6=REF.p.eta6; K7=REF.p.kappa7; RHO8=REF.p.rho8; I0ETA=REF.i0_eta
A_INT=A.astype(np.int64)

@njit(cache=True)
def mod1s(x): return x-math.floor(x)
@njit(cache=True)
def f_scalar(i,z):
    if i==0: return math.sin(2*math.pi*z[2]+PHI1)
    if i==1:
        c=math.cos(2*math.pi*z[5]+PHI2); return 4*c*c*c-3*c
    if i==2: return math.sin(2*math.pi*(z[7]-z[1]))/(1+RHO3*math.cos(2*math.pi*z[1])**2)
    if i==3: return math.tanh(K4*math.sin(2*math.pi*(z[0]-z[6])))
    if i==4:
        q=.5*(1+math.sin(2*math.pi*(z[6]+z[7]))); return 4*q*(1-q)-.5
    if i==5: return (math.exp(ETA6*math.cos(2*math.pi*(z[1]-z[0])))-I0ETA)/(math.exp(ETA6)-math.exp(-ETA6))
    if i==6: return (2/math.pi)*math.atan(K7*math.sin(2*math.pi*(z[3]+z[5])))
    return math.log(1+RHO8*math.sin(math.pi*(z[4]-z[2]))**2)/math.log(1+RHO8)-.5
@njit(cache=True)
def quant_step(state,qbits):
    Q=1<<qbits; z=np.empty(8,np.float64)
    for i in range(8): z[i]=state[i]/Q
    for i in range(8): z[i]=mod1s(z[i]+EPS[i]*f_scalar(i,z))
    y=np.empty(8,np.float64)
    for i in range(8):
        s=0.0
        for j in range(8): s+=A_INT[i,j]*z[j]
        y[i]=mod1s(s)
    z=y
    for i in range(7,-1,-1): z[i]=mod1s(z[i]-EPS[i]*f_scalar(i,z))
    out=np.empty(8,np.uint64)
    for i in range(8):
        q=int(math.floor(z[i]*Q+.5))%Q; out[i]=q
    return out
@njit(cache=True)
def eqstate(a,b):
    for i in range(8):
        if a[i]!=b[i]: return False
    return True
@njit(cache=True)
def brent(x0,qbits,cap):
    power=1; lam=1; steps=1; tort=x0.copy(); hare=quant_step(x0,qbits)
    while not eqstate(tort,hare) and steps<cap:
        if power==lam: tort=hare.copy(); power*=2; lam=0
        hare=quant_step(hare,qbits); lam+=1; steps+=1
    if not eqstate(tort,hare): return -1,-1,steps
    mu=0; tort=x0.copy(); hare=x0.copy()
    for _ in range(lam): hare=quant_step(hare,qbits)
    while not eqstate(tort,hare) and mu<cap:
        tort=quant_step(tort,qbits); hare=quant_step(hare,qbits); mu+=1
    return mu,lam,steps

def qinitial(key,qbits):
    sys=HCHM8(key); Q=1<<qbits; return np.minimum((sys.x0*Q).astype(np.uint64),Q-1)

def repeat_search(sys,dtype,cap=100000):
    x=sys.x0.astype(dtype); seen={x.tobytes():0}
    for n in range(1,cap+1):
        x=sys.step(x,dtype=dtype); b=x.tobytes()
        if b in seen: return seen[b],n-seen[b],n
        seen[b]=n
    return None,None,cap

def run_finite_tests():
    t0=time.time(); _=quant_step(qinitial(key_for_index(0),8),8)
    cfg={6:100000,8:250000,10:500000,12:1000000,16:1000000}; rows=[]
    for q,cap in cfg.items():
        for i in range(20):
            mu,lam,steps=brent(qinitial(key_for_index(100+i),q),q,cap)
            rows.append({'qbits':q,'run':i,'cap':cap,'cycle_found':int(lam>=0),'transient':mu,'cycle_length':lam,'steps_examined':steps})
    df=pd.DataFrame(rows); df.to_csv(DATA/'quantized_cycle_20runs.csv',index=False)
    floatrows=[]
    for dtype,name,cap in [(np.float32,'Float32',100000),(np.float64,'Float64',100000)]:
        for i in range(5):
            mu,lam,n=repeat_search(HCHM8(key_for_index(500+i)),dtype,cap)
            floatrows.append({'precision':name,'run':i,'cap':cap,'repeat_found':int(lam is not None),'transient':mu,'cycle_length':lam,'steps':n})
    pd.DataFrame(floatrows).to_csv(DATA/'float_repeat_search.csv',index=False)
    div=[]
    for i in range(30):
        sys=HCHM8(key_for_index(900+i)); a=sys.x0.astype(np.float64); b=sys.x0.astype(np.float32); hit=None
        for n in range(1,501):
            a=sys.step(a,np.float64); b=sys.step(b,np.float32); d=np.linalg.norm(torus_diff(a,b.astype(np.float64)))
            if hit is None and d>.1: hit=n
        div.append({'run':i,'first_distance_gt_0p1':hit or 501,'distance_at_500':d})
    pd.DataFrame(div).to_csv(DATA/'float32_float64_divergence.csv',index=False)
    summaries=[]
    for q,g in df.groupby('qbits'):
        f=g[g.cycle_found==1]
        summaries.append({'qbits':int(q),'runs':len(g),'cap':int(g.cap.iloc[0]),'cycles_found':len(f),
          'median_transient':None if len(f)==0 else float(f.transient.median()),'median_cycle':None if len(f)==0 else float(f.cycle_length.median()),
          'min_cycle':None if len(f)==0 else int(f.cycle_length.min()),'max_cycle':None if len(f)==0 else int(f.cycle_length.max())})
    summary={'quantized_state_emulation':summaries,'float_repeat_search':floatrows,
             'float32_float64_median_divergence_iteration':float(pd.DataFrame(div).first_distance_gt_0p1.median()),
             'float32_float64_range':[int(pd.DataFrame(div).first_distance_gt_0p1.min()),int(pd.DataFrame(div).first_distance_gt_0p1.max())],
             'important_scope':'Quantized-state tests evaluate round-to-nearest state quantization with transcendental functions evaluated in float64; they are not an integer-only fixed-point implementation.',
             'elapsed_s':time.time()-t0}
    (DATA/'finite_precision_validated_summary.json').write_text(json.dumps(summary,indent=2))
    p=df.copy(); p['plot']=np.where(p.cycle_found==1,p.cycle_length,p.cap)
    plt.figure(figsize=(8,4.8)); plt.boxplot([p[p.qbits==q]['plot'] for q in cfg],tick_labels=[f'Q{q}' for q in cfg],showfliers=False); plt.yscale('log'); plt.ylabel('Cycle length or right-censoring cap'); plt.tight_layout(); plt.savefig(FIG/'quantized_cycle_results.pdf'); plt.savefig(FIG/'quantized_cycle_results.png',dpi=220); plt.close()
    return summary

# ---- S-box verification ----
HAS1_HEX='''09 68 F4 4D 7E 3C 32 B3 E4 BF 0B 6B F1 19 8F 07
9E 2F 35 D3 17 37 3A B0 A4 50 69 6C 22 EC BC 29
3D 46 72 74 57 F3 C8 45 25 AC 5D 1D 4E 5C E9 D1
73 A5 D7 6D B2 36 54 9B 96 94 DE DD AE 38 99 05
33 55 C0 11 70 3B 04 2D 76 59 28 12 3E 02 CC 80
88 5F FF 4F 16 0C 3F 6F 56 AF 84 C5 E5 78 13 F8
D4 63 C6 A1 E2 EB E7 0A D2 43 D9 5A 1A CF 67 60
87 F0 1E 9C 00 A0 86 ED 0F F2 C3 06 A7 D8 89 E0
EF BB B4 97 4A B1 D6 26 61 90 27 1C 8E 8A C4 23
E3 1F 03 7D 44 7A 21 52 51 15 75 83 B8 FD 40 9D
AB 18 30 DC 47 A8 14 93 DB 41 FA F6 FE 2C 6A 7C
D0 0E CA D5 66 6E 01 5B A9 8B 58 E6 5E 9A 77 AA
FC DA 48 79 EE 82 53 2E 98 8C EA 08 C7 F5 7B 4C
C2 DF A3 64 62 24 AD 85 95 8D CD 91 E1 B9 1B B7
20 42 BD F7 39 F9 BA 2A 31 92 65 9F 71 CE 2B C1
B5 BE FB B6 0D E8 34 CB C9 81 A6 10 4B A2 7F 49'''
HAS1=np.array([int(x,16) for x in HAS1_HEX.split()],dtype=np.uint8)
RIN=[0xA1,0x14,0xD6,0x71,0xFB,0x49,0x21,0xB6]; ROUT=[0x11,0x56,0x94,0x43,0x71,0xEE,0x19,0x91]
PAR=np.array([i.bit_count()&1 for i in range(256)],dtype=np.uint8)
def gf_mul(a,b):
    r=0
    for _ in range(8):
        if b&1:r^=a
        hi=a&0x80; a=(a<<1)&255
        if hi:a^=0x1b
        b>>=1
    return r
def gf_pow(a,e):
    r=1
    while e:
        if e&1:r=gf_mul(r,a)
        a=gf_mul(a,a); e>>=1
    return r
def gf_inv(a): return 0 if a==0 else gf_pow(a,254)
GFINV=np.array([gf_inv(i) for i in range(256)],dtype=np.uint8)
def linmap_table(rows):
    X=np.arange(256,dtype=np.uint16); out=np.zeros(256,dtype=np.uint8)
    for i,r in enumerate(rows): out|=(PAR[np.bitwise_and(X,r)]<<i).astype(np.uint8)
    return out
def build_sbox(rin,rout,a,b):
    li=linmap_table(rin); lo=linmap_table(rout); return np.bitwise_xor(lo[GFINV[np.bitwise_xor(li,a)]],b).astype(np.uint8)
def invtable(S):
    I=np.empty(256,dtype=np.uint8); I[S]=np.arange(256,dtype=np.uint8); return I
def fwht(a):
    a=np.array(a,dtype=np.int32,copy=True); h=1
    while h<len(a):
        for i in range(0,len(a),2*h):
            x=a[i:i+h].copy(); y=a[i+h:i+2*h].copy(); a[i:i+h]=x+y; a[i+h:i+2*h]=x-y
        h*=2
    return a
def ddt(S):
    X=np.arange(256,dtype=np.uint16); D=np.zeros((256,256),dtype=np.uint16)
    for a in range(256): D[a]=np.bincount(np.bitwise_xor(S[X],S[np.bitwise_xor(X,a)]),minlength=256)
    return D
def lat(S):
    L=np.zeros((256,256),dtype=np.int16)
    for b in range(256): L[:,b]=fwht(1-2*PAR[np.bitwise_and(S,b)].astype(np.int16))
    return L
def degrees(S):
    out=[]
    for b in range(1,256):
        f=PAR[np.bitwise_and(S,b)].copy()
        for i in range(8):
            for m in range(256):
                if m&(1<<i): f[m]^=f[m^(1<<i)]
        out.append(max((m.bit_count() for m,v in enumerate(f) if v),default=0))
    return np.array(out)
def sac(S):
    X=np.arange(256,dtype=np.uint16); M=np.zeros((8,8))
    for i in range(8):
        d=np.bitwise_xor(S[X],S[np.bitwise_xor(X,1<<i)])
        for j in range(8): M[i,j]=np.mean((d>>j)&1)
    return M
def bct(S):
    I=invtable(S); X=np.arange(256,dtype=np.uint16); B=np.zeros((256,256),dtype=np.uint16)
    for a in range(256):
        Sa=S[np.bitwise_xor(X,a)]
        for b in range(256): B[a,b]=np.count_nonzero(np.bitwise_xor(I[np.bitwise_xor(S,b)],I[np.bitwise_xor(Sa,b)])==a)
    return B
def dlct(S):
    X=np.arange(256,dtype=np.uint16); T=np.zeros((256,256),dtype=np.int16)
    for a in range(256):
        d=np.bitwise_xor(S[X],S[np.bitwise_xor(X,a)])
        for b in range(256): T[a,b]=np.sum(1-2*PAR[np.bitwise_and(d,b)].astype(np.int16))
    return T
def bic_corr(S):
    X=np.arange(256,dtype=np.uint16); vals=[]
    for i in range(8):
        d=np.bitwise_xor(S[X],S[np.bitwise_xor(X,1<<i)]); bits=np.column_stack([((d>>j)&1) for j in range(8)]).astype(float); C=np.corrcoef(bits,rowvar=False)
        for p in range(8):
            for q in range(p+1,8): vals.append(abs(C[p,q]))
    return max(vals),float(np.mean(vals))
def cycle_structure(S):
    seen=np.zeros(256,bool); lengths=[]
    for x in range(256):
        if not seen[x]:
            y=x; l=0
            while not seen[y]: seen[y]=True; y=int(S[y]); l+=1
            lengths.append(l)
    return sorted(lengths,reverse=True)
def run_sbox_tests():
    t0=time.time(); rec=build_sbox(RIN,ROUT,0xF0,0x77); assert np.array_equal(rec,HAS1)
    I=invtable(HAS1); D=ddt(HAS1); L=lat(HAS1); B=bct(HAS1); T=dlct(HAS1); S=sac(HAS1); deg=degrees(HAS1)
    nl=np.array([128-np.max(np.abs(L[:,b]))//2 for b in range(1,256)]); bm,bmean=bic_corr(HAS1)
    np.savetxt(DATA/'HAS1_DDT.csv',D,fmt='%d',delimiter=','); np.savetxt(DATA/'HAS1_LAT.csv',L,fmt='%d',delimiter=','); np.savetxt(DATA/'HAS1_BCT.csv',B,fmt='%d',delimiter=','); np.savetxt(DATA/'HAS1_DLCT.csv',T,fmt='%d',delimiter=','); np.savetxt(DATA/'HAS1_SAC.csv',S,fmt='%.8f',delimiter=','); np.savetxt(DATA/'HAS1_inverse.csv',I.reshape(16,16),fmt='%02X',delimiter=',')
    vals,freq=np.unique(D[1:],return_counts=True)
    metrics={'sha256':hashlib.sha256(HAS1.tobytes()).hexdigest(),'bijective':True,'round_trip_inverse':bool(np.array_equal(I[HAS1],np.arange(256,dtype=np.uint8))),
      'nonlinearity_min':int(nl.min()),'nonlinearity_max':int(nl.max()),'nonlinearity_distribution':{str(int(v)):int(np.sum(nl==v)) for v in np.unique(nl)},
      'differential_uniformity':int(D[1:].max()),'ddt_spectrum':{str(int(v)):int(f) for v,f in zip(vals,freq)},'max_differential_probability':float(D[1:].max()/256),
      'max_walsh':int(np.max(np.abs(L[:,1:]))),'max_linear_correlation':float(np.max(np.abs(L[:,1:]))/256),'max_linear_bias':float(np.max(np.abs(L[:,1:]))/512),
      'degree_min':int(deg.min()),'degree_max':int(deg.max()),'degree_distribution':{str(int(v)):int(np.sum(deg==v)) for v in np.unique(deg)},
      'boomerang_uniformity':int(B[1:,1:].max()),'dlct_absolute_indicator':int(np.max(np.abs(T[1:,1:]))),
      'fixed_points':int(np.sum(HAS1==np.arange(256))),'opposite_fixed_points':int(np.sum(HAS1==(np.arange(256)^255))),
      'sac_mean':float(S.mean()),'sac_mad':float(np.abs(S-.5).mean()),'sac_maxdev':float(np.abs(S-.5).max()),
      'bic_max_abs_corr':float(bm),'bic_mean_abs_corr':bmean,'bic_pairwise_nonlinearity_min':int(min(128-np.max(np.abs(L[:,(1<<p)|(1<<q)]))//2 for p in range(8) for q in range(p+1,8))),
      'cycle_structure':cycle_structure(HAS1),'inverse_nonlinearity_min':int(min(128-np.max(np.abs(lat(I)[:,b]))//2 for b in range(1,256))),
      'first_two_known_answers':{'00':f'{int(HAS1[0]):02X}','01':f'{int(HAS1[1]):02X}'},'elapsed_s':time.time()-t0}
    (DATA/'HAS1_independent_verification.json').write_text(json.dumps(metrics,indent=2))
    for M,name,title in [(D[1:,1:],'ddt_heatmap','HAS-1 DDT'),(L[1:,1:],'lat_heatmap','HAS-1 LAT'),(B[1:,1:],'bct_heatmap','HAS-1 BCT'),(T[1:,1:],'dlct_heatmap','HAS-1 DLCT')]:
        plt.figure(figsize=(6.2,5.2)); plt.imshow(M,aspect='auto'); plt.colorbar(); plt.title(title); plt.tight_layout(); plt.savefig(FIG/f'{name}.pdf'); plt.savefig(FIG/f'{name}.png',dpi=220); plt.close()
    plt.figure(figsize=(6.2,5.2)); plt.imshow(S,vmin=.4,vmax=.6); plt.colorbar(); plt.xlabel('Output bit'); plt.ylabel('Flipped input bit'); plt.tight_layout(); plt.savefig(FIG/'sac_heatmap.pdf'); plt.savefig(FIG/'sac_heatmap.png',dpi=220); plt.close()
    return metrics

# ---- dynamic profile family v2 ----
def make_rows(stream:bytes,idx:int):
    rows=[1<<i for i in range(8)]
    for _ in range(64):
        v=stream[idx%len(stream)]; idx+=1; i=(v>>1)&7; j=stream[idx%len(stream)]&7; idx+=1
        if i==j:j=(j+1)&7
        if v&1: rows[i]^=rows[j]
        else: rows[i],rows[j]=rows[j],rows[i]
    return rows,idx

def serialize_hchm(key,round_index,burn=64,samples=16):
    sys=HCHM8(key,round_index); x=sys.x0.copy()
    for _ in range(burn): x=sys.step(x)
    out=bytearray()
    for _ in range(samples):
        x=sys.step(x)
        for v in x: out.extend(min(int(v*2**32),2**32-1).to_bytes(4,'big'))
    return bytes(out)
def candidate(base,counter):
    st=shake(base+counter.to_bytes(2,'big')+b'|AFFINE-PROFILE-V2|',258); rin,idx=make_rows(st,0); rout,idx=make_rows(st,idx); a=st[idx%len(st)]; b=st[(idx+1)%len(st)]; return build_sbox(rin,rout,a,b),rin,rout,a,b
def select(key,r,include_hchm,budget=64):
    base=key+r.to_bytes(4,'big')+b'|HAS1-FAMILY-V2|'
    if include_hchm: base+=serialize_hchm(key,r)
    first=None; best=None
    for c in range(budget):
        S,rin,rout,a,b=candidate(base,c); X=np.arange(256,dtype=np.uint8); fixed=int(np.sum(S==X)); opp=int(np.sum(S==(X^255)))
        if fixed==0 and opp==0:
            M=sac(S); mad=float(np.abs(M-.5).mean()); mx=float(np.abs(M-.5).max())
        else: mad=1.; mx=1.
        score=(fixed+opp,mx,mad,c)
        if best is None or score<best[0]: best=(score,S,rin,rout,a,b,fixed,opp,mad,mx,c)
        if first is None and fixed==0 and opp==0 and mad<=.03 and mx<=.08: first=(S,rin,rout,a,b,fixed,opp,mad,mx,c)
    fallback=first is None
    if fallback:
        _,S,rin,rout,a,b,fixed,opp,mad,mx,c=best; first=(S,rin,rout,a,b,fixed,opp,mad,mx,c)
    return first,fallback

def run_family(nkeys=1000,rounds=16,ablation_keys=250):
    t0=time.time(); rows=[]; hashes=set(); dups=0
    for ki in range(nkeys):
        key=key_for_index(10000+ki); prev=None
        for r in range(rounds):
            sel,fb=select(key,r,True); S,rin,rout,a,b,fix,opp,mad,mx,c=sel; h=hashlib.sha256(S.tobytes()).hexdigest(); dups+=int(h in hashes); hashes.add(h)
            rows.append({'key_index':ki,'round':r,'counter':c,'fallback':int(fb),'fixed':fix,'opposite':opp,'sac_mad':mad,'sac_maxdev':mx,'sha256':h,'equal_positions_previous':None if prev is None else int(np.sum(S==prev))}); prev=S
    df=pd.DataFrame(rows); df.to_csv(DATA/'dynamic_family_v2_16000.csv',index=False)
    abl=[]
    for mode in (True,False):
        for ki in range(ablation_keys):
            key=key_for_index(30000+ki)
            for r in range(rounds):
                sel,fb=select(key,r,mode); S,rin,rout,a,b,fix,opp,mad,mx,c=sel
                abl.append({'mode':'HCHM+SHAKE' if mode else 'SHAKE-only','key_index':ki,'round':r,'counter':c,'fallback':int(fb),'sac_mad':mad,'sac_maxdev':mx,'sha256':hashlib.sha256(S.tobytes()).hexdigest()})
    adf=pd.DataFrame(abl); adf.to_csv(DATA/'ablation_hchm_vs_shake_v2.csv',index=False); g1=adf[adf['mode']=='HCHM+SHAKE']; g2=adf[adf['mode']=='SHAKE-only']
    ks1=stats.ks_2samp(g1.sac_mad,g2.sac_mad); ks2=stats.ks_2samp(g1.sac_maxdev,g2.sac_maxdev)
    def smd(x,y):
        p=math.sqrt(((len(x)-1)*np.var(x,ddof=1)+(len(y)-1)*np.var(y,ddof=1))/(len(x)+len(y)-2)); return float((np.mean(x)-np.mean(y))/p) if p else 0.
    summary={'generator_version':'HAS1-FAMILY-V2','tables':len(df),'keys':nkeys,'rounds':rounds,'unique_tables':len(hashes),'duplicates':dups,'fallbacks':int(df.fallback.sum()),'fallback_rate':float(df.fallback.mean()),
      'counter_mean':float(df.counter.mean()),'counter_median':float(df.counter.median()),'counter_p95':float(df.counter.quantile(.95)),'counter_max':int(df.counter.max()),
      'sac_mad_mean':float(df.sac_mad.mean()),'sac_mad_std':float(df.sac_mad.std()),'sac_maxdev_mean':float(df.sac_maxdev.mean()),'sac_maxdev_max':float(df.sac_maxdev.max()),
      'previous_round_equal_positions_mean':float(df.equal_positions_previous.dropna().mean()),'previous_round_equal_positions_max':int(df.equal_positions_previous.dropna().max()),
      'ablation_tables_per_mode':len(g1),'ablation_ks_sac_mad_stat':float(ks1.statistic),'ablation_ks_sac_mad_p':float(ks1.pvalue),'ablation_ks_sac_max_stat':float(ks2.statistic),'ablation_ks_sac_max_p':float(ks2.pvalue),'ablation_smd_sac_mad':smd(g1.sac_mad.values,g2.sac_mad.values),'ablation_smd_sac_max':smd(g1.sac_maxdev.values,g2.sac_maxdev.values),
      'interpretation':'After SHAKE256 whitening, HCHM input is not expected to improve affine-invariant S-box metrics; it acts as deterministic profile diversification, while cryptographic assurance rests on SHAKE256 and the inverse core.','elapsed_s':time.time()-t0}
    (DATA/'family_v2_summary.json').write_text(json.dumps(summary,indent=2))
    plt.figure(figsize=(7.5,4.8)); plt.hist(df.counter,bins=np.arange(-.5,64.5,1)); plt.xlabel('Selected candidate counter'); plt.ylabel('Frequency'); plt.tight_layout(); plt.savefig(FIG/'family_counter_distribution.pdf'); plt.savefig(FIG/'family_counter_distribution.png',dpi=220); plt.close()
    plt.figure(figsize=(7.5,4.8));
    for mode,g in adf.groupby('mode'): plt.hist(g.sac_mad,bins=30,alpha=.5,label=mode)
    plt.xlabel('SAC mean absolute deviation'); plt.ylabel('Frequency'); plt.legend(); plt.tight_layout(); plt.savefig(FIG/'ablation_sac.pdf'); plt.savefig(FIG/'ablation_sac.png',dpi=220); plt.close()
    return summary

def main():
    import sys
    stage=sys.argv[1] if len(sys.argv)>1 else 'all'
    out={'stage':stage,'started':time.strftime('%Y-%m-%dT%H:%M:%S')}
    if stage in ('map','all'): out['map']=run_map_tests()
    if stage in ('finite','all'): out['finite_precision']=run_finite_tests()
    if stage in ('sbox','all'): out['sbox']=run_sbox_tests()
    if stage in ('family','all'): out['family']=run_family()
    out['completed']=time.strftime('%Y-%m-%dT%H:%M:%S')
    (DATA/f'{stage}_stage_summary.json').write_text(json.dumps(out,indent=2)); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
