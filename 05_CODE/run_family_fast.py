import sys,time,json,hashlib,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
sys.path.insert(0,str(Path(__file__).resolve().parent))
import run_validated_experiments as m
DATA=Path(__file__).resolve().parents[1]/'data'; FIG=Path(__file__).resolve().parents[1]/'figures'

def direct_family(nkeys=256,rounds=16):
    rows=[]; hashes=set(); duplicates=0; t0=time.time()
    for mode in (True,False):
        for ki in range(nkeys):
            key=m.key_for_index(40000+ki); prev=None
            for r in range(rounds):
                base=key+r.to_bytes(4,'big')+b'|HAS1-FAMILY-V2|'
                if mode: base+=m.serialize_hchm(key,r)
                S,rin,rout,a,b=m.candidate(base,0)
                h=hashlib.sha256(S.tobytes()).hexdigest(); duplicates+=int((mode,h) in hashes); hashes.add((mode,h))
                X=np.arange(256,dtype=np.uint8); fixed=int(np.sum(S==X)); opp=int(np.sum(S==(X^255)))
                equal=None if prev is None else int(np.sum(S==prev)); prev=S
                # SAC on stratified 1/8 subset only
                if (ki*rounds+r)%8==0:
                    M=m.sac(S); mad=float(np.abs(M-.5).mean()); mx=float(np.abs(M-.5).max())
                else: mad=np.nan; mx=np.nan
                rows.append({'mode':'HCHM+SHAKE' if mode else 'SHAKE-only','key_index':ki,'round':r,'counter':0,'fixed':fixed,'opposite':opp,'sac_mad':mad,'sac_maxdev':mx,'sha256':h,'equal_positions_previous':equal})
    df=pd.DataFrame(rows); df.to_csv(DATA/'dynamic_family_direct_8192.csv',index=False)
    # selected fixed-budget subset 32 keys x 4 rounds per mode
    selrows=[]
    for mode in (True,False):
        for ki in range(32):
            key=m.key_for_index(50000+ki)
            for r in range(4):
                sel,fb=m.select(key,r,mode,budget=64); S,rin,rout,a,b,fix,opp,mad,mx,c=sel
                selrows.append({'mode':'HCHM+SHAKE' if mode else 'SHAKE-only','key_index':ki,'round':r,'counter':c,'fallback':int(fb),'fixed':fix,'opposite':opp,'sac_mad':mad,'sac_maxdev':mx,'sha256':hashlib.sha256(S.tobytes()).hexdigest()})
    sdf=pd.DataFrame(selrows); sdf.to_csv(DATA/'fixed_budget_subset_256.csv',index=False)
    # stats
    h=df[df['mode']=='HCHM+SHAKE']; s=df[df['mode']=='SHAKE-only']; hs=h.dropna(subset=['sac_mad']); ss=s.dropna(subset=['sac_mad'])
    ks=stats.ks_2samp(hs.sac_mad,ss.sac_mad)
    summary={'direct_profiles_per_mode':len(h),'total_direct_profiles':len(df),'unique_hashes_total':df.sha256.nunique(),'duplicates_within_mode':duplicates,
      'HCHM_fixed_mean':float(h.fixed.mean()),'SHAKE_fixed_mean':float(s.fixed.mean()),'HCHM_opposite_mean':float(h.opposite.mean()),'SHAKE_opposite_mean':float(s.opposite.mean()),
      'HCHM_consecutive_equal_positions_mean':float(h.equal_positions_previous.dropna().mean()),'SHAKE_consecutive_equal_positions_mean':float(s.equal_positions_previous.dropna().mean()),
      'sac_subset_per_mode':len(hs),'HCHM_sac_mad_mean':float(hs.sac_mad.mean()),'SHAKE_sac_mad_mean':float(ss.sac_mad.mean()),'KS_sac_mad_stat':float(ks.statistic),'KS_sac_mad_p':float(ks.pvalue),
      'fixed_budget_subset_per_mode':int(len(sdf)/2),'HCHM_fallback_rate':float(sdf[sdf['mode']=='HCHM+SHAKE'].fallback.mean()),'SHAKE_fallback_rate':float(sdf[sdf['mode']=='SHAKE-only'].fallback.mean()),
      'HCHM_counter_median':float(sdf[sdf['mode']=='HCHM+SHAKE'].counter.median()),'SHAKE_counter_median':float(sdf[sdf['mode']=='SHAKE-only'].counter.median()),
      'interpretation':'All generated S-boxes are affine-equivalent to the inverse function, so nonlinearity 112, DU 4, degree 7 and BU 6 are invariant. The corpus evaluates profile diversity and implementation-selection behavior, not a new equivalence class.',
      'elapsed_s':time.time()-t0}
    (DATA/'family_fast_summary.json').write_text(json.dumps(summary,indent=2))
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7.2,4.5));
    for mode,g in df.dropna(subset=['sac_mad']).groupby('mode'): plt.hist(g.sac_mad,bins=25,alpha=.5,label=mode)
    plt.legend(); plt.xlabel('SAC mean absolute deviation'); plt.ylabel('Frequency'); plt.tight_layout(); plt.savefig(FIG/'family_ablation_sac.pdf'); plt.savefig(FIG/'family_ablation_sac.png',dpi=220); plt.close()
    plt.figure(figsize=(7.2,4.5));
    for mode,g in sdf.groupby('mode'): plt.hist(g.counter,bins=np.arange(-.5,64.5,1),alpha=.5,label=mode)
    plt.legend(); plt.xlabel('Selected counter in 64-candidate budget'); plt.ylabel('Frequency'); plt.tight_layout(); plt.savefig(FIG/'fixed_budget_counter_subset.pdf'); plt.savefig(FIG/'fixed_budget_counter_subset.png',dpi=220); plt.close()
    print(json.dumps(summary,indent=2))
if __name__=='__main__': direct_family()
