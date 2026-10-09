import sys,os,json,pathlib
S=pathlib.Path(__file__).resolve().parent;R=pathlib.Path(__file__).resolve().parents[4]
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']:os.environ[k]=str(S/k)
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']: (S/k).mkdir(exist_ok=True)
sys.path.insert(0,str(R/'pipeline'))
import env,lcsearch as L,numpy as np,pandas as pd
from scipy.ndimage import median_filter
from astropy.timeseries import BoxLeastSquares
print('PID',os.getpid(),'port none',flush=True);out={}
for tic in [370612437,252721481]:
 raw,_,_=L.reject_bad_sectors(L.load_star(tic));t,f,sec=raw['time'],raw['flux'],raw['sector'];flat=np.full(len(t),np.nan);ee=np.zeros(len(t))
 for s in np.unique(sec):
  ix=np.where(sec==s)[0];cuts=np.r_[0,np.where(np.diff(t[ix])>.5)[0]+1,len(ix)]
  for a,b in zip(cuts[:-1],cuts[1:]):
   j=ix[a:b];size=int(1./np.median(np.diff(t[j])))//2*2+1;trend=median_filter(f[j],size=size,mode='nearest');flat[j]=f[j]/trend
   flat[j[(t[j]-t[j[0]]<.5)|(t[j[-1]]-t[j]<.5)]]=np.nan
  sd=L.robust_std(flat[ix][np.isfinite(flat[ix])]);ee[ix]=sd
 keep=np.isfinite(flat)&(abs(flat-1)<5*ee);c={'time':t[keep],'flux':flat[keep],'flux_err':ee[keep],'sector':sec[keep],'exptime':raw['exptime'][keep]};c=L.bin_lc(c)
 tg=pd.read_csv(R/'results/targets_ranked.csv').set_index('TIC').loc[tic];rho=tg.Mass/tg.Rad**3
 periods=L.period_grid(np.ptp(c['time']),rho);bls=BoxLeastSquares(c['time'],c['flux'],c['flux_err']);res=bls.power(periods,np.array([.5,.75,1,1.5,2,3,4])/24,objective='snr');sde,_=L.sde(res.power,periods)
 k=int(np.argmax(sde));kr=int(np.argmax(res.power));out[str(tic)]={'best_SDE_period':float(res.period[k]),'SDE':float(sde[k]),'SNR':float(res.depth_snr[k]),'best_raw_SNR_period':float(res.period[kr]),'n_periods':len(periods)}
 (S/'search.json').write_text(json.dumps(out,indent=2));print(tic,json.dumps(out[str(tic)]),flush=True)
