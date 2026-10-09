import os,sys,json, pathlib, requests
S=pathlib.Path(__file__).resolve().parent; R=pathlib.Path(__file__).resolve().parents[4]
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']:os.environ[k]=str(S/k)
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']: (S/k).mkdir(exist_ok=True)
sys.path.insert(0,str(R/'pipeline'))
import env,vet as V, numpy as np,pandas as pd,lightkurve as lk
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import importlib
M=importlib.import_module('10_pixels')
print('PID',os.getpid(),'port none',flush=True);out={}
for tic in [370612437,252721481]:
 tg=pd.read_csv(R/'results/targets_ranked.csv').set_index('TIC').loc[tic]
 query=f'SELECT source_id,ra,dec,pmra,pmdec FROM gaiadr3.gaia_source WHERE source_id={int(tg.GAIA)}'
 cache=S/f'gaia_{tic}.csv'
 if not cache.exists():
  r=requests.get('https://gea.esac.esa.int/tap-server/tap/sync',params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query},timeout=40)
  r.raise_for_status(); cache.write_text(r.text)
 me=pd.read_csv(cache).iloc[0]
 P,t0,dur=V.ephemeris(tic); saved=json.loads((V.OUT/f'{tic}_pixels.json').read_text());secs=[];alls=[]
 for sector in saved['sectors']:
  sec=sector['sector']; hits=sorted((R/'cache/tesscut').rglob(f'tess-s{sec:04d}-*_{tg.ra:.6f}_{tg.dec:.6f}_*.fits'))
  tpf=lk.read(str(hits[0]),quality_bitmask='hard');t=np.array(tpf.time.value);cube=np.array(tpf.flux.value,dtype=float)
  ok=np.isfinite(t)&np.isfinite(cube).all(axis=(1,2));t=t[ok];cube=cube[ok];direct=np.median(cube,axis=0);faint=direct<np.percentile(direct,30);cube-=np.median(cube[:,faint],axis=1)[:,None,None]
  tc=SkyCoord(*M.at_epoch(me,np.median(t)+2457000),unit='deg');tx,ty=map(float,tpf.wcs.world_to_pixel(tc));yy,xx=np.mgrid[:15,:15];near=np.hypot(xx-tx,yy-ty)<4.5
  imgs=M.diff_images(t,cube,P,t0,dur); nulls=np.array([M.diff_images(t,cube,P,t0+o*P,dur).mean(0) for o in V.NULL_OFFSETS]);noise=1.4826*np.median(abs(nulls-np.median(nulls,axis=0)),axis=0)
  s={'yy':yy,'xx':xx,'near':near,'noise':noise,'sigma':sector['psf_sigma_px'],'target_px':(tx,ty)};secs.append(s);alls.append(np.concatenate([imgs.mean(0)[None],nulls]));print(tic,sec,'ready',flush=True)
  del cube,tpf
 scores=[]
 for j in range(41):
  diffs=[a[j] for a in alls];scores.append(M.chi2_none(secs,diffs)-M.chi2_at(secs,diffs,[s['target_px'] for s in secs])[0])
 out[str(tic)]={'observed':scores[0],'null_scores':scores[1:],'n_null_ge_real':sum(np.array(scores[1:])>=scores[0]),'empirical_p_add_one':(1+sum(np.array(scores[1:])>=scores[0]))/41}
 (S/'pixel_null.json').write_text(json.dumps(out,indent=2,default=lambda x:x.item()));print(tic,json.dumps(out[str(tic)],default=lambda x:x.item()),flush=True)
