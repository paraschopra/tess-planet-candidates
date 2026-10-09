import os,sys,json,pathlib
S=pathlib.Path(__file__).resolve().parent;R=pathlib.Path(__file__).resolve().parents[4]
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']:os.environ[k]=str(S/k)
for k in ['MPLCONFIGDIR','XDG_CACHE_HOME','XDG_CONFIG_HOME','ASTROPY_CACHE_DIR','NUMBA_CACHE_DIR']: (S/k).mkdir(exist_ok=True)
sys.path.insert(0,str(R/'pipeline'))
import env,vet as V,lcsearch as L,numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from lightkurve.utils import TessQualityFlags
print('PID',os.getpid(),'port none',flush=True)
def score(t,y,ex,P,t0,dur,pc=None):
 ep=np.rint((t-t0)/P).astype(int);num=den=0
 for e in np.unique(ep):
  x=t-t0-e*P;m=abs(x)<.3;xx=x[m]
  if sum(xx<-dur)<4 or sum(xx>dur)<4:continue
  dt=ex/86400;box=np.maximum(0,np.minimum(xx+dt/2,dur/2)-np.maximum(xx-dt/2,-dur/2))/dt
  if box.sum()<.8:continue
  A=np.column_stack([np.ones(len(xx)),xx]) if pc is None else np.column_stack([np.ones(len(xx)),xx,pc[m]-np.mean(pc[m],axis=0)])
  residual=box-A@np.linalg.lstsq(A,box,rcond=None)[0];num+=residual@y[m];den+=residual@residual
 return num/den if den>0 else np.nan
out={}
for tic in [370612437,252721481]:
 f=json.loads((V.OUT/f'{tic}_fit.json').read_text())['stellar'];P,t0,dur=f['P']['median'],f['t0']['median'],f['t14_hr']['median']/24
 offsets=np.r_[0,V.NULL_OFFSETS];allval=[];allcorr=[];pers=[]
 for path in sorted((R/'cache/lightkurve').rglob(f'*{tic:016d}-s*_lc.fits')):
  with fits.open(path) as h:
   d=h[1].data;sec=h[0].header['SECTOR'];exp=1800 if sec<27 else 200
   m=((d['QUALITY']&TessQualityFlags.DEFAULT_BITMASK)==0)&np.isfinite(d['TIME'])&np.isfinite(d['MOM_CENTR1'])&np.isfinite(d['MOM_CENTR2'])&np.isfinite(d['POS_CORR1'])&np.isfinite(d['POS_CORR2'])
   t=np.array(d['TIME'][m],float);c=np.column_stack([d['MOM_CENTR1'][m],d['MOM_CENTR2'][m]]);pc=np.column_stack([d['POS_CORR1'][m],d['POS_CORR2'][m]])
   w=WCS(h[2].header);z=w.pixel_to_world(5.,5.)
   jac=np.array([[a.to_value(u.mas) for a in z.spherical_offsets_to(w.pixel_to_world(6.,5.))],[a.to_value(u.mas) for a in z.spherical_offsets_to(w.pixel_to_world(5.,6.))]]).T
  arr=np.array([[score(t,c[:,k],exp,P,t0+o*P,dur) for k in [0,1]] for o in offsets]); cor=np.array([[score(t,c[:,k],exp,P,t0+o*P,dur,pc=pc) for k in [0,1]] for o in offsets])
  sky=arr@jac.T;skyc=cor@jac.T;allval.append(sky);allcorr.append(skyc)
  pers.append({'sector':sec,'raw_millipixel':(arr[0]*1000).tolist(),'raw_sky_mas':sky[0].tolist()});print(tic,sec,'ready',flush=True)
 rr={'per_sector':pers}
 for name,v in [('raw',allval),('pointing_regressed',allcorr)]:
  a=np.array(v);np.save(S/f'centroid_{tic}_{name}.npy',a);err=1.4826*np.median(abs(a[:,1:]-np.median(a[:,1:],axis=1)[:,None]),axis=1);weights=1/err**2;combined=np.sum(a*weights[:,None,:],axis=0)/np.sum(weights,axis=0);sd=1.4826*np.median(abs(combined[1:]-np.median(combined[1:],axis=0)),axis=0)
  rr[name]={'sky_shift_mas':combined[0].tolist(),'control_MAD_mas':sd.tolist(),'z':(combined[0]/sd).tolist()}
 out[str(tic)]=rr;(S/'centroid.json').write_text(json.dumps(out,indent=2));print(tic,json.dumps(rr),flush=True)
