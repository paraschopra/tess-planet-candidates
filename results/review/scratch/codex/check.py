import os, sys, json, time
from pathlib import Path
SCR=Path(__file__).resolve().parent; ROOT=Path(__file__).resolve().parents[4]
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']: os.environ[k]=str(SCR/k)
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']: (SCR/k).mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'pipeline'))
import env, lcsearch as L, vet as V
import numpy as np, pandas as pd
from scipy.stats import chi2
from astropy.io import fits
import lightkurve as lk
import importlib
DL=importlib.import_module('02_download')
rng=np.random.default_rng(90713)
def local(raw,P,t0,dur):
    # Linear baseline + exposure-integrated box, fit independently around each predicted event.
    t,f=raw['time'],raw['flux']; epoch=np.rint((t-t0)/P).astype(int); rows=[]
    for s in np.unique(raw['sector']):
        m=raw['sector']==s; noise=L.ptp_noise(f[m])
        for ep in np.unique(epoch[m]):
            x=t-t0-ep*P; take=m&(abs(x)<0.3); xx=x[take]; yy=f[take]
            exp=raw['exptime'][take]/86400
            if sum(xx<-dur)<4 or sum(xx>dur)<4: continue
            box=np.maximum(0,np.minimum(xx+exp/2,dur/2)-np.maximum(xx-exp/2,-dur/2))/exp
            if box.sum()<0.8: continue
            A=np.column_stack([np.ones(len(xx)),xx,-box]); cov=np.linalg.inv(A.T@A)*noise**2
            b=np.linalg.lstsq(A,yy,rcond=None)[0]; rows.append([int(s),int(ep),float(b[2]),float(cov[2,2])])
    a=np.array(rows); w=1/a[:,3]; d=np.sum(w*a[:,2])/w.sum(); err=1/np.sqrt(w.sum())
    boot=[]
    for _ in range(3000):
        ix=rng.integers(0,len(a),len(a)); ww=w[ix]; boot.append(np.sum(ww*a[ix,2])/ww.sum())
    return {'n_events':len(a),'depth_ppm':d*1e6,'white_error_ppm':err*1e6,'event_bootstrap_error_ppm':float(np.std(boot)*1e6),'positive_events':int(sum(a[:,2]>0)), 'season_depths':[{'sectors':list(map(int,g)),'depth':float(np.sum(w[np.isin(a[:,0],g)]*a[np.isin(a[:,0],g),2])/sum(w[np.isin(a[:,0],g)])*1e6),'error':float(1/np.sqrt(sum(w[np.isin(a[:,0],g)]))*1e6)} for g in V.seasons(raw)]}
if __name__ == '__main__':
    out={'pid':os.getpid(),'port':None,'scratch':str(SCR)}; print(json.dumps(out),flush=True)
    for tic in [370612437,252721481]:
        raw,_,_=L.reject_bad_sectors(L.load_star(tic)); P,t0,dur=V.ephemeris(tic); fit=json.loads((V.OUT/f'{tic}_fit.json').read_text())
        rr={'local_phase2':local(raw,P,t0,dur),'local_fitted':local(raw,fit['stellar']['P']['median'],fit['stellar']['t0']['median'],fit['stellar']['t14_hr']['median']/24)}
        print(tic,'LOCAL',json.dumps(rr),flush=True)
        # Rerun matched processing at original and two interleaved phase grids.
        offsets=np.concatenate([V.NULL_OFFSETS, np.linspace(.095,.405,40),np.linspace(.595,.905,40)])
        depths=[]
        for i,o in enumerate(offsets):
            c=V.masked_clean(raw,P,t0+o*P,dur)
            depths.append(V.box_depth(c,abs(V.phase_of(c['time'],P,t0+o*P))<dur/2)*1e6)
            if i%40==39: print(tic,'controls',i+1,flush=True)
        ds=np.array(depths); sig=lambda v:float(1.4826*np.median(abs(v-np.median(v))))
        rr['controls']={'old40_MAD':sig(ds[:40]),'new80_MAD':sig(ds[40:]),'all120_MAD':sig(ds),'all120_SD':float(np.std(ds,ddof=1)), 'depths':depths}
        rr['time_range']=list(map(float,[min(raw['time']),max(raw['time'])]))
        qparts=[]
        for path in sorted((ROOT/'cache/lightkurve').rglob(f'*{tic:016d}*_llc.fits')):
            q=lk.read(str(path),quality_bitmask='default'); t,f,e,col=DL.flux_arrays(q,'QLP'); sec=int(q.meta['SECTOR']); exp=1800 if sec<27 else 600 if sec<56 else 200
            qparts.append({'time':t,'flux':f,'flux_err':e,'sector':np.full(len(t),sec),'exptime':np.full(len(t),exp)})
        qraw={k:np.concatenate([x[k] for x in qparts]) for k in qparts[0]}; order=np.argsort(qraw['time']);qraw={k:v[order] for k,v in qraw.items()};qraw,_,_=L.reject_bad_sectors(qraw)
        d,e,null=V.depth_emp_matched(qraw,P,t0,dur);rr['qlp_matched']={'depth_ppm':d*1e6,'error_ppm':e*1e6,'snr':d/e}
        out[str(tic)]=rr;(SCR/'checks.json').write_text(json.dumps(out,indent=2));print(tic,'DONE',json.dumps({k:v for k,v in rr.items() if k!='controls'}),flush=True)
