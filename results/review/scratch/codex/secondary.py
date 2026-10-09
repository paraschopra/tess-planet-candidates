from check import *
out={}
for tic in [370612437,252721481]:
    raw,_,_=L.reject_bad_sectors(L.load_star(tic))
    f=json.loads((V.OUT/f'{tic}_fit.json').read_text())['stellar']
    P,t0,d=f['P']['median'],f['t0']['median'],f['t14_hr']['median']/24
    vals=[]
    for pp,tt in [(P,t0+P/2),(2*P,t0),(2*P,t0+P)]:
        r=local(raw,pp,tt,d)
        vals.append([r['depth_ppm'],r['event_bootstrap_error_ppm']])
    out[str(tic)]={'secondary_odd_even':vals,'odd_even_z':float(abs(vals[1][0]-vals[2][0])/np.hypot(vals[1][1],vals[2][1]))}
print(json.dumps(out,indent=2))
(SCR/'secondary.json').write_text(json.dumps(out,indent=2))
