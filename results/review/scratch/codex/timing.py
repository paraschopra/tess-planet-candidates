import os,json,pathlib,sys
S=pathlib.Path(__file__).resolve().parent
os.environ['XDG_CACHE_HOME']=str(S/'XDG_CACHE_HOME');os.environ['XDG_CONFIG_HOME']=str(S/'XDG_CONFIG_HOME')
for k in ['XDG_CACHE_HOME','XDG_CONFIG_HOME']: (S/k).mkdir(exist_ok=True)
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[4]/'pipeline'))
import env
import pandas as pd,numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord,EarthLocation
from astropy.utils import iers
iers.conf.auto_download=False
import astropy.units as u
R=pathlib.Path(__file__).resolve().parents[4];tg=pd.read_csv(R/'results/targets_ranked.csv').set_index('TIC');out={}
for tic in [370612437,252721481]:
 d=pd.read_csv(R/f'results/phase3/{tic}_predicted_transits.csv'); f=json.loads((R/f'results/phase3/{tic}_fit.json').read_text())['stellar'];correct=f['t0']['median']+2457000+d.epoch*f['P']['median'];err_sec=(d.mid_BJD_TDB-correct)*86400
 sky=SkyCoord(tg.loc[tic,'ra'],tg.loc[tic,'dec'],unit='deg');loc=EarthLocation.from_geocentric(0,0,0,unit='m');bary=Time(d.mid_BJD_TDB.values,format='jd',scale='tdb',location=loc);guess=bary
 for i in range(4):guess=bary-guess.light_travel_time(sky,kind='barycentric')
 original=Time(d.mid_UTC.tolist(),scale='utc');shift=(original-guess.utc).to_value(u.min)
 out[str(tic)]={'n':len(d),'max_BJD_resid_seconds':float(abs(err_sec).max()),'labelled_UTC_minus_geocentric_UTC_min':[float(shift.min()),float(shift.max())],'first_row':d.iloc[0].to_dict(),'correct_geocentric_UTC_first':guess.utc.iso[0]}
print(json.dumps(out,indent=2));(S/'timing.json').write_text(json.dumps(out,indent=2))
