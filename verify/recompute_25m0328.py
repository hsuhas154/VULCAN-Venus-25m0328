#!/usr/bin/env python3
"""Recompute every quantity quoted in Experiments 7, 8 and 9 from committed outputs.
Emits a registry (JSON) plus a readable table. Read-only."""
import os, sys, pickle, json, numpy as np
from scipy import interpolate

if not os.path.isdir("output") or not os.path.isdir("atm"):
    sys.exit("ERROR: run this from the repository root, where output/ and atm/ live.")

R="output/"
def L(p):
    with open(R+p,"rb") as f: return pickle.load(f)

A1  = L("Experiment2_A1_UUV-nominal-25m0328.vul")
A0  = L("Experiment1_A0_UUV-control-25m0328.vul")
E7A = L("Experiment7A_UUV-upward18-25m0328.vul")
E7B = L("Experiment7B_UUV-downward18-25m0328.vul")
B8  = L("Experiment8_baseline_A1-25m0328.vul")
E8A = L("Experiment8A_ClS2-25m0328.vul")
E8B = L("Experiment8B_OSSO_S2O2-25m0328.vul")
E8C = L("Experiment8C_S3-25m0328.vul")
E9  = L("Experiment9_ClS2-selfconsistent-25m0328.vul")

z    = np.asarray(A1["atm"]["zmco"])/1e5
dzcm = np.asarray(A1["atm"]["dz"]); dzkm = dzcm/1e5
bins = np.asarray(A1["variable"]["bins"])
quad = np.where(bins<240.,0.1,2.0)
band = (z>=69)&(z<=93)
AREA = np.pi*(0.98e-4)**2/4.0   # pi D^2 / 4 with D = D_1, as op.py forms it
TAU_REQ = 1.094968472081
REG={}

def put(k,v):
    REG[k]=float(v); return v

def radiative(case, ref, tag):
    d=np.asarray(case["variable"]["aflux"])-np.asarray(ref["variable"]["aflux"])
    f0=np.asarray(ref["variable"]["aflux"])
    norm=np.sqrt(np.sum(f0**2*quad[None,:]))
    ij=np.unravel_index(np.argmax(np.abs(d)),d.shape)
    l2=np.sqrt(np.sum(d**2*quad[None,:]))
    l1z=np.sum(np.abs(d)*quad[None,:],axis=1); w=l1z*dzkm; p=w/w.sum()
    l2z=np.sqrt(np.sum(d**2*quad[None,:],axis=1))
    l2b=np.sqrt(np.sum(f0**2*quad[None,:],axis=1))
    frac=np.where(l2b>0,l2z/l2b,0.0); j=int(np.argmax(frac))
    cq=np.cumsum(p); qs=[float(z[np.searchsorted(cq,q)]) for q in (.25,.5,.75,.90)]
    put(f"{tag}.peak",d[ij]); put(f"{tag}.peak_z",z[ij[0]]); put(f"{tag}.peak_lam",bins[ij[1]])
    put(f"{tag}.L2pct",100*l2/norm); put(f"{tag}.meanalt",np.sum(z*p))
    put(f"{tag}.inband",100*p[band].sum()); put(f"{tag}.p5968",100*p[(z>=59)&(z<=68)].sum())
    put(f"{tag}.above93",100*p[z>93].sum()); put(f"{tag}.below59",100*p[z<59].sum())
    put(f"{tag}.maxfrac",100*frac[j]); put(f"{tag}.maxfrac_z",z[j])
    for i,q in enumerate(qs): put(f"{tag}.q{[25,50,75,90][i]}",q)
    return d

d7A=radiative(E7A,A1,"e7a"); d7B=radiative(E7B,A1,"e7b")
d6 =radiative(A1,A0,"e6ctl")
d8A=radiative(E8A,B8,"e8a"); d8B=radiative(E8B,B8,"e8b"); d8C=radiative(E8C,B8,"e8c")
d9 =radiative(E9 ,B8,"e9")

# --- Exp 7 derived claims
put("e7.ratio_7A_over_7B", REG["e7a.L2pct"]/REG["e7b.L2pct"])
put("e7.meanalt_shift", abs(REG["e7a.meanalt"]-REG["e7b.meanalt"]))

# --- profiles
def prof(path,tag):
    d=np.loadtxt(path,skiprows=1); zk=d[:,0]/1e5 if d[:,0].max()>1e4 else d[:,0]
    n=d[:,1]; t=n.sum()
    put(f"{tag}.levelsum",t); put(f"{tag}.centroid",(zk*n).sum()/t)
    put(f"{tag}.peak_z",zk[np.argmax(n)]); put(f"{tag}.peak_N",n.max())
    put(f"{tag}.inband",100*n[(zk>=69)&(zk<=93)].sum()/t)
    put(f"{tag}.p5869",100*n[(zk>=58)&(zk<69)].sum()/t)
    N=np.interp(z,zk,n); put(f"{tag}.colhi69",np.sum(N[z>=69]*dzcm[z>=69]))
    put(f"{tag}.col_dzw",np.sum(N*dzcm))
    return zk,n
prof("atm/mode1+2.txt","pnom")
prof("atm/mode1+2_Experiment9_ClS2profile_25m0328.txt","p9")
put("e9.centroid_shift", REG["pnom.centroid"]-REG["p9.centroid"])

# --- retained tau
grid=bins[(bins>=300)&(bins<=698)]
def tau69(spec,pf,tag,slin):
    uv=np.loadtxt(spec,skiprows=1)
    q=interpolate.interp1d(uv[:,0],uv[:,1],kind="slinear" if slin else "linear",
                           bounds_error=False,fill_value=0.0)(grid)
    m=(grid>=300)&(grid<=620)
    d=np.loadtxt(pf,skiprows=1); zk=d[:,0]/1e5 if d[:,0].max()>1e4 else d[:,0]
    N=np.interp(z,zk,d[:,1]); hi=z>=69
    put(f"{tag}.bandint",np.sum(q*2.0))
    put(f"{tag}.peakQ",q.max()); put(f"{tag}.peakQ_lam",grid[np.argmax(q)])
    sup=grid[q>0]
    put(f"{tag}.sup_lo",sup.min()); put(f"{tag}.sup_hi",sup.max())
    return put(f"{tag}.tau69",np.sum(q[m]*2.0)*AREA*float(np.sum(N[hi]*dzcm[hi])))
tau69("atm/UV_absorber.txt","atm/mode1+2.txt","spec_emp",True)
tau69("atm/UV_absorber_Experiment8A_ClS2_25m0328.txt","atm/mode1+2.txt","spec_8A",False)
tau69("atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt","atm/mode1+2.txt","spec_8B",False)
tau69("atm/UV_absorber_Experiment8C_S3_25m0328.txt","atm/mode1+2.txt","spec_8C",False)
tau69("atm/UV_absorber_Experiment8A_ClS2_25m0328.txt",
      "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt","spec_9",False)

# --- opacity screen (exact reproduction of the committed script)
var,atm=A1["variable"],A1["atm"]
y=np.asarray(var["y"]); ymix=np.asarray(var["ymix"]); species=list(var["species"])
above=z>=58.0; bandw=(bins>=300)&(bins<=620)
CAND=["SCl2","ClS2","S2Cl2","SCl","S3","S4","S2","S2O","S2O2","SO3","Cl2","ClO","OSCl"]
screen={}
for nm in CAND:
    if nm not in species or nm not in var["cross"]: continue
    sig=np.asarray(var["cross"][nm]); col=float(np.sum(y[:,species.index(nm)][above]*dzcm[above]))
    t=np.where(bandw,sig*col,0.0); j=int(np.argmax(t))
    if t[j]<=0: continue
    screen[nm]=(bins[j],sig[j],col,t[j],TAU_REQ/t[j])
    put(f"scr.{nm}.lam",bins[j]); put(f"scr.{nm}.sigma",sig[j])
    put(f"scr.{nm}.col",col); put(f"scr.{nm}.tau",t[j]); put(f"scr.{nm}.short",TAU_REQ/t[j])
sig=np.asarray(var["cross"]["ClS2"]); colC=float(np.sum(y[:,species.index("ClS2")][above]*dzcm[above]))
j=int(np.argmax(np.where(bandw,sig*colC,0.0))); needed=TAU_REQ/sig[j]
put("cls2.needed_col",needed)
hcl=float(np.sum(y[:,species.index("HCl")][above]*dzcm[above]))
put("cls2.frac_HCl",needed/hcl)
layer=(z>=58)&(z<=70); path=float(np.sum(dzcm[layer]))
meanM=float(np.average(np.asarray(atm["M"])[layer],weights=dzcm[layer]))
put("cls2.req_mix",needed/path/meanM)
put("hcl.mix60",ymix[np.argmin(np.abs(z-60)),species.index("HCl")])
put("cls2.sup_lo",bins[sig>0].min()); put("cls2.sup_hi",bins[sig>0].max())
put("cls2.sup_peak",bins[np.argmax(sig)])

# --- spectral correlation with the empirical curve, and 320-400 opacity share
uv=np.loadtxt("atm/UV_absorber.txt",skiprows=1)
qe=interpolate.interp1d(uv[:,0],uv[:,1],kind="slinear",bounds_error=False,fill_value=0.0)(grid)
for tag,f in [("8A","atm/UV_absorber_Experiment8A_ClS2_25m0328.txt"),
              ("8B","atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt"),
              ("8C","atm/UV_absorber_Experiment8C_S3_25m0328.txt")]:
    u=np.loadtxt(f,skiprows=1)
    qq=interpolate.interp1d(u[:,0],u[:,1],kind="linear",bounds_error=False,fill_value=0.0)(grid)
    m=(grid>=300)&(grid<=620)
    put(f"corr.{tag}",np.corrcoef(qe[m],qq[m])[0,1])
    sh=(grid>=320)&(grid<=400)
    put(f"share320400.{tag}",100*np.sum(qq[sh]*2.)/np.sum(qq[m]*2.))
m=(grid>=300)&(grid<=620); sh=(grid>=320)&(grid<=400)
put("share320400.emp",100*np.sum(qe[sh]*2.)/np.sum(qe[m]*2.))

# --- photolysis
PAR=["S2O2","SCl2","ClS2","S3","S2Cl2","SCl","S2O","S4","SO3","H2SO4","SO2","SO"]
def jmax(case,ref,p,mode):
    k=(p,0)
    if k not in ref["variable"]["J_sp"]: return None
    A=np.asarray(ref["variable"]["J_sp"][k]); B=np.asarray(case["variable"]["J_sp"][k])
    if mode=="dJ": a,b=A,B
    else:
        i=list(ref["variable"]["species"]).index(p)
        a=A*np.asarray(ref["variable"]["y"])[:,i]; b=B*np.asarray(case["variable"]["y"])[:,i]
    d=b-a; j=int(np.argmax(np.where(band,np.abs(d),0.0)))
    return (100*d[j]/a[j] if a[j] else float("nan")), z[j]
for tag,case,ref in [("e8a",E8A,B8),("e8b",E8B,B8),("e8c",E8C,B8),("e9",E9,B8),("e7b",E7B,A1)]:
    for p in PAR:
        r=jmax(case,ref,p,"dJ")
        if r: put(f"dJ.{tag}.{p}",r[0]); put(f"dJz.{tag}.{p}",r[1])
for p in PAR:
    r=jmax(E9,B8,p,"dJy")
    if r: put(f"dJy.e9.{p}",r[0]); put(f"dJyz.e9.{p}",r[1])

# --- reservoirs
RES=["SO2","SO","SO3","H2SO4","H2SO4_l"]
def resv(case,ref,tag):
    sp=list(ref["variable"]["species"])
    yb=np.asarray(ref["variable"]["y"]); yc=np.asarray(case["variable"]["y"])
    nb=np.asarray(ref["atm"]["n_0"]); nc=np.asarray(case["atm"]["n_0"])
    for s in RES:
        if s not in sp: continue
        i=sp.index(s); mb=yb[:,i]/nb; mc=yc[:,i]/nc
        ok=(mb>1e-20)&band; rel=np.zeros_like(mb)
        rel[ok]=np.abs(mc[ok]-mb[ok])/mb[ok]
        cb=np.sum(yb[:,i]*dzcm); cc=np.sum(yc[:,i]*dzcm)
        put(f"res.{tag}.{s}.inband",100*rel.max())
        put(f"res.{tag}.{s}.col",100*(cc-cb)/cb)
for tag,c in [("e8a",E8A),("e8b",E8B),("e8c",E8C),("e9",E9)]: resv(c,B8,tag)

# --- convergence
for tag,s in [("b8",B8),("e8a",E8A),("e8b",E8B),("e8c",E8C),("e9",E9)]:
    put(f"conv.{tag}.longdy",s["variable"]["longdy"])

# --- Exp 9 identity metrics
s=d9+d6
put("id.max_sum",np.abs(s).max())
put("id.max_sum_pct",100*np.abs(s).max()/np.abs(d6).max())
l2=lambda x: np.sqrt(np.sum(x**2*quad[None,:]))
put("id.l2_ratio_pct",100*l2(s)/l2(d6))
bm=band[:,None]*np.ones_like(d6,dtype=bool)
put("id.l2_band_pct",100*l2(np.where(bm,s,0))/l2(np.where(bm,d6,0)))
with np.errstate(divide='ignore', invalid='ignore'):
    r=np.where(np.abs(d6)>0.01*np.abs(d6).max(),-d9/d6,np.nan)
put("id.ratio_mean",np.nanmean(r)); put("id.ratio_sd",np.nanstd(r))
put("id.ratio_min",np.nanmin(r)); put("id.ratio_max",np.nanmax(r))
put("id.ncells",np.sum(~np.isnan(r)))

# --- channel ratios e9/e7b
rs=[]
for p in PAR:
    a=REG.get(f"dJ.e9.{p}"); b=REG.get(f"dJ.e7b.{p}")
    if a is None or b is None or abs(b)<=0.1: continue
    rs.append(a/b); put(f"chr.{p}",a/b)
rs=np.array(rs)
put("chr.n",len(rs)); put("chr.mean",rs.mean()); put("chr.sd",rs.std())
put("chr.spread",rs.max()-rs.min())
i81=int(np.argmin(np.abs(z-81)))
put("chr.flux_ratio",np.sum(np.abs(d9[i81])*quad)/np.sum(np.abs(d7B[i81])*quad))

# --- ClS2 mixing ratio at model levels
mr=np.asarray(A1["variable"]["y"])[:,species.index("ClS2")]/np.asarray(A1["atm"]["n_0"])
for zz in [0,55,57,59,61,63,65,67,69]:
    k=int(np.argmin(np.abs(z-zz))); put(f"cls2mix.{zz}",mr[k])
put("cls2mix.peak",mr.max()); put("cls2mix.col58",colC)
# cross-check ymix against y/n_0
put("xcheck.ymix_vs_y_over_n0", float(np.max(np.abs(
    ymix[:,species.index("ClS2")] - mr))))


# ============================================================================
# CORRECTED CONVENTIONS.  op.py line 62 interpolates the particle profile onto
# zco[1:], the upper interface of each layer, with slinear, and the radiative
# transfer weights it by dz.  Every absorber-column quantity uses that and only
# that.  Bands are named constants; none is written inline.
# ============================================================================
BAND_RESP   = (69.0, 93.0)
BAND_LOWER  = (35.0, 59.0)
BAND_CLOUD  = (58.0, 68.999)
GRID_NM     = (300.0, 698.0)
SCREEN_NM   = (300.0, 620.0)

zco = np.asarray(A1["atm"]["zco"])/1e5
zu  = zco[1:]
assert len(zu) == len(dzcm), "interface grid and dz length disagree"

def Nprof(path):
    d = np.loadtxt(path, skiprows=1)
    return interpolate.interp1d(d[:,0], d[:,1], kind="slinear",
                                bounds_error=False, fill_value=0.0)(zu)

def profile_stats(path, tag):
    N = Nprof(path)
    col = float(np.sum(N*dzcm))
    put(f"{tag}.levelsum", N.sum())
    put(f"{tag}.column",   col)
    put(f"{tag}.centroid", np.sum(zu*N*dzcm)/col)
    put(f"{tag}.peak_z",   zu[np.argmax(N)])
    put(f"{tag}.peak_N",   N.max())
    for nm,(lo,hi) in [("inband",BAND_RESP),("lower",BAND_LOWER),("cloud",BAND_CLOUD)]:
        m=(zu>=lo)&(zu<=hi); put(f"{tag}.{nm}", 100*np.sum(N[m]*dzcm[m])/col)
    m=(zu>=0)&(zu<=30); put(f"{tag}.below30", 100*np.sum(N[m]*dzcm[m])/col)
    hi=zu>=BAND_RESP[0]; put(f"{tag}.colhi69", np.sum(N[hi]*dzcm[hi]))
    put(f"{tag}.lowest_share", 100*N[0]/N.sum())
    return N

profile_stats("atm/mode1+2.txt","pnom")
profile_stats("atm/mode1+2_Experiment9_ClS2profile_25m0328.txt","p9")
profile_stats("atm/mode1+2_Experiment7A_shift18_25m0328.txt","p7a")
profile_stats("atm/mode1+2_Experiment7B_shiftDown18_25m0328.txt","p7b")
put("p9.col_deficit_pct", 100*(1-REG["p9.column"]/REG["pnom.column"]))

def Qgrid(path, slin):
    u = np.loadtxt(path, skiprows=1)
    return interpolate.interp1d(u[:,0], u[:,1], kind="slinear" if slin else "linear",
                                bounds_error=False, fill_value=0.0)(grid)

def tau_above(spec, pf, tag, slin):
    q = Qgrid(spec, slin); m = (grid>=SCREEN_NM[0]) & (grid<=SCREEN_NM[1])
    N = Nprof(pf); hi = zu >= BAND_RESP[0]
    put(f"{tag}.bandint", np.sum(q*2.0))
    put(f"{tag}.peakQ", q.max()); put(f"{tag}.peakQ_lam", grid[np.argmax(q)])
    sup = grid[q>0]; put(f"{tag}.sup_lo", sup.min()); put(f"{tag}.sup_hi", sup.max())
    return put(f"{tag}.tau69", float(np.sum(q[m]*2.0))*AREA*float(np.sum(N[hi]*dzcm[hi])))

tau_above("atm/UV_absorber.txt","atm/mode1+2.txt","spec_emp",True)
for t,f in [("spec_8A","atm/UV_absorber_Experiment8A_ClS2_25m0328.txt"),
            ("spec_8B","atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt"),
            ("spec_8C","atm/UV_absorber_Experiment8C_S3_25m0328.txt")]:
    tau_above(f,"atm/mode1+2.txt",t,False)
tau_above("atm/UV_absorber_Experiment8A_ClS2_25m0328.txt",
          "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt","spec_9",False)
put("tau.reduction", REG["spec_emp.tau69"]/REG["spec_9.tau69"])
put("tau.spec_change_pct", 100*(REG["spec_8A.tau69"]-REG["spec_emp.tau69"])/REG["spec_emp.tau69"])

# scale factors, from the model cross sections as the generator does
I0 = float(np.sum(Qgrid("atm/UV_absorber.txt",True)*2.0))
gm = (bins>=GRID_NM[0]) & (bins<=GRID_NM[1])
for t,sp in [("8A","ClS2"),("8B","S2O2"),("8C","S3")]:
    sig = np.asarray(A1["variable"]["cross"][sp])[gm]
    put(f"k.{t}", I0/float(np.sum(sig*2.0)))

# spectral similarity and band share, both over the FULL 300 to 698 nm grid
qe = Qgrid("atm/UV_absorber.txt",True)
sh = (grid>=320)&(grid<=400)
put("share320400.emp", 100*np.sum(qe[sh]*2.)/np.sum(qe*2.))
put("share320400.emp_620", 100*np.sum(qe[sh]*2.)/np.sum(qe[(grid>=300)&(grid<=620)]*2.))
for t,f in [("8A","atm/UV_absorber_Experiment8A_ClS2_25m0328.txt"),
            ("8B","atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt"),
            ("8C","atm/UV_absorber_Experiment8C_S3_25m0328.txt")]:
    q = Qgrid(f,False)
    put(f"corr.{t}", np.corrcoef(qe,q)[0,1])
    put(f"share320400.{t}", 100*np.sum(q[sh]*2.)/np.sum(q*2.))

# constants that appear as literals
put("const.tau_req", TAU_REQ); put("const.area", AREA)
put("const.D_UV", 0.98e-4)
for tag,s_ in [("a1",A1),("e7a",E7A),("e7b",E7B),("b8",B8),("e8a",E8A),
               ("e8b",E8B),("e8c",E8C),("e9",E9)]:
    for k in ("longdy","longdydt","dt"):
        if k in s_["variable"]: put(f"conv.{tag}.{k}", s_["variable"][k])



# ---- Experiment 7 tables -----------------------------------------------
SPEC_BANDS=[(300,320),(320,340),(340,360),(360,380),(380,420),(420,500),(500,600),(600,698)]
for tag,dd in [("e7a",d7A),("e7b",d7B)]:
    w=np.abs(dd)*quad[None,:]*dzkm[:,None]; tot=w.sum()
    for lo,hi in SPEC_BANDS:
        m=(bins>=lo)&(bins<=hi) if hi==698 else (bins>=lo)&(bins<hi)
        put(f"specband.{tag}.{lo}_{hi}", 100*w[:,m].sum()/tot)
for p in PAR:
    k=(p,0)
    if k not in A1["variable"]["J_sp"]: continue
    J=np.asarray(A1["variable"]["J_sp"][k])
    put(f"Jpeak.full.{p}", J.max())
    put(f"Jpeak.inband.{p}", J[band].max())
    for tag,c in [("e7a",E7A),("e7b",E7B)]:
        B_=np.asarray(c["variable"]["J_sp"][k]); dd=B_-J
        j=int(np.argmax(np.where(band,np.abs(dd),0.0)))
        put(f"dJ.{tag}.{p}", 100*dd[j]/J[j]); put(f"dJz.{tag}.{p}", z[j])
for tag,c in [("e7a",E7A),("e7b",E7B)]: resv(c,A1,tag)
put("corr.emp",1.0)
put("id.inv_l2_ratio", 100.0/REG["id.l2_ratio_pct"])
put("e9.tau_spec_change_pct", 100*(REG["spec_8A.tau69"]-REG["spec_emp.tau69"])/REG["spec_emp.tau69"])


# ---- Experiment 10 -----------------------------------------------------
E10 = L("Experiment10_S3-on-ClS2profile-25m0328.vul")
d10 = radiative(E10, B8, "e10")
for k in ("longdy","longdydt","dt"): put(f"conv.e10.{k}", E10["variable"][k])

tau_above("atm/UV_absorber.txt","atm/mode1+2_Experiment7A_shift18_25m0328.txt","spec_7A",True)
tau_above("atm/UV_absorber.txt","atm/mode1+2_Experiment7B_shiftDown18_25m0328.txt","spec_7B",True)
tau_above("atm/UV_absorber_Experiment8C_S3_25m0328.txt",
          "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt","spec_10",False)
put("tau.7A_over_A1", REG["spec_7A.tau69"]/REG["spec_emp.tau69"])
put("tau.9_minus_10", abs(REG["spec_9.tau69"]-REG["spec_10.tau69"]))

_L2 = lambda x: float(np.sqrt(np.sum(x**2*quad[None,:])))
_diff = d10 - d9
put("e10.fielddiff_pct", 100*_L2(_diff)/_L2(d9))
for nm,(lo,hi_) in [("inband",(69,93)),("c5968",(59,68)),("above93",(94,112)),("below59",(0,58))]:
    mz=(z>=lo)&(z<=hi_)
    put(f"e10.fd_{nm}", 100*_L2(np.where(mz[:,None],_diff,0))/_L2(np.where(mz[:,None],d9,0)))
    put(f"e9.norm_{nm}", _L2(np.where(mz[:,None],d9,0)))

PAR8=["S2O2","SCl2","ClS2","S3","S2Cl2","SCl","S2O","S4"]
for p in PAR8:
    r=jmax(E10,B8,p,"dJ")
    if r: put(f"dJ.e10.{p}", r[0]); put(f"dJz.e10.{p}", r[1])
def _ratio(tagA,tagB,name):
    r=np.array([REG[f"dJ.{tagA}.{p}"]/REG[f"dJ.{tagB}.{p}"] for p in PAR8])
    put(f"{name}.mean", r.mean()); put(f"{name}.spread", r.max()-r.min())
    put(f"{name}.sd", r.std()); return r
_ratio("e10","e7b","chr8_10_7b"); _ratio("e9","e7b","chr8_9_7b")
_ratio("e10","e9","chr8_10_9");  _ratio("e7a","e7b","chr8_7a_7b")
_i81=int(np.argmin(np.abs(z-81))); _bw=(bins>=300)&(bins<=620)
_d7b=np.asarray(E7B["variable"]["aflux"])-np.asarray(A1["variable"]["aflux"])
put("chr8.flux_ratio_10_7b", float(np.sum(d10[_i81]*quad*_bw)/np.sum(_d7b[_i81]*quad*_bw)))
put("chr8.flux_ratio_9_7b",  float(np.sum(d9[_i81]*quad*_bw)/np.sum(_d7b[_i81]*quad*_bw)))
put("chr8.agree_10_7b", abs(REG["chr8_10_7b.mean"]-REG["chr8.flux_ratio_10_7b"]))
put("chr8.agree_9_7b",  abs(REG["chr8_9_7b.mean"]-REG["chr8.flux_ratio_9_7b"]))

d8a=np.asarray(E8A["variable"]["aflux"])-np.asarray(B8["variable"]["aflux"])
d8c=np.asarray(E8C["variable"]["aflux"])-np.asarray(B8["variable"]["aflux"])
put("twoby2.nominal_ratio", _L2(d8c)/_L2(d8a))
put("twoby2.cls2_ratio",    _L2(d10)/_L2(d9))
put("twoby2.nominal_fielddiff_pct", 100*_L2(d8c-d8a)/_L2(d8a))
resv(E10,B8,"e10")

REG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "registry.json")
json.dump(REG, open(REG_PATH, "w"), indent=1)
print(f"registry written to {REG_PATH}")
print(f"registry entries: {len(REG)}")
for k in sorted(REG):
    print(f"  {k:32s} {REG[k]:.6g}")
