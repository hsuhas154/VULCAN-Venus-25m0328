#!/usr/bin/env python3
"""Recompute every quantity quoted in Experiments 1 through 13 from committed outputs.
Emits a registry (JSON) plus a readable table. Read-only: no VULCAN rerun.

TWO PHOTOLYSIS CONVENTIONS, WHICH MUST NEVER BE MIXED
-----------------------------------------------------
    dJ       the change in the photolysis RATE COEFFICIENT. This is the purely
             radiative response, and it peaks where the flux change peaks.
    d(J y)   the change in the photolysis RATE itself. This folds in the parent
             species abundance and therefore peaks LOWER in the atmosphere,
             because abundance rises going down while flux falls.

Experiment 6 reported d(J y). Experiments 7 through 10 reported dJ. A number
carried from one convention into a sentence written for the other will look like
an altitude disagreement of several kilometres and is not one. This note
supersedes the standalone verify_all_25m0328.py, which recorded it first.

FOUR COMPARISON CONVENTIONS, RECOVERED FROM THE DOCUMENTS THEMSELVES
--------------------------------------------------------------------
 1. A relative change divides by max(|reference|, |case|), not by the reference.
    Experiment 2's actinic-flux section prints this formula explicitly; the
    abundance tables of Experiments 1, 2 and 3 use it without saying so.
 2. A mean over a saved field is unweighted, over every element of the array.
    The "bounded relative" means instead average only over elements whose
    denominator is non-zero.
 3. A band integral carries VULCAN's endpoint half-weighting at the band edges.
    Omitting it inflates every band by a few per cent.
 4. Where a tendency is integrated over altitude, the layer thickness must be in
    cm for the result to be a column rate in cm^-2 s^-1. Experiment 7's table was
    printed with the thickness in km, which leaves it 1e5 low; both forms are
    registered, as .col and .col_km.
"""
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


# ---- Experiment 11 -----------------------------------------------------
E11A = L("Experiment11A_SSCl2-25m0328.vul"); E11B = L("Experiment11B_ClSSCl-25m0328.vul")
radiative(E11A, B8, "e11a"); radiative(E11B, B8, "e11b")
for t,s_ in [("e11a",E11A),("e11b",E11B)]:
    for k in ("longdy","longdydt","dt"): put(f"conv.{t}.{k}", s_["variable"][k])
tau_above("atm/UV_absorber_Experiment11A_SSCl2_25m0328.txt","atm/mode1+2.txt","spec_11A",False)
tau_above("atm/UV_absorber_Experiment11B_ClSSCl_25m0328.txt","atm/mode1+2.txt","spec_11B",False)
_m620=(grid>=300)&(grid<=620)
put("spec_emp.int620", float(np.sum(Qgrid("atm/UV_absorber.txt",True)[_m620]*2.0)))
for t,f in [("11A","atm/UV_absorber_Experiment11A_SSCl2_25m0328.txt"),
            ("11B","atm/UV_absorber_Experiment11B_ClSSCl_25m0328.txt")]:
    q=Qgrid(f,False)
    put(f"corr.{t}", np.corrcoef(qe,q)[0,1])
    put(f"share320400.{t}", 100*np.sum(q[sh]*2.)/np.sum(q*2.))
for t,c in [("e11a",E11A),("e11b",E11B)]:
    for p in PAR:
        r=jmax(c,B8,p,"dJ")
        if r: put(f"dJ.{t}.{p}", r[0]); put(f"dJz.{t}.{p}", r[1])
    resv(c,B8,t)

# isomer opacity screen with the computed cross sections
_col_s2cl2 = float(np.sum(np.asarray(A1["variable"]["y"])[:,
    list(A1["variable"]["species"]).index("S2Cl2")][z>=58.0]*dzcm[z>=58.0]))
put("iso.col_s2cl2", _col_s2cl2)
_PUB={"SSCl2":(3.37e-17,264.0),"ClSSCl":(2.63e-17,240.0)}
for tag,f in [("SSCl2","csec-sscl2.dat"),("ClSSCl","csec-clsscl.dat")]:
    d_=np.loadtxt("Experiment11/data/"+f); o=np.argsort(d_[:,0])
    lam,sig=d_[o,0],d_[o,1]; mm=(lam>=SCREEN_NM[0])&(lam<=SCREEN_NM[1])
    gp,gl=sig.max(),lam[int(np.argmax(sig))]
    bp,bl=sig[mm].max(),lam[mm][int(np.argmax(sig[mm]))]
    put(f"iso.{tag}.global_sigma",gp); put(f"iso.{tag}.global_lam",gl)
    put(f"iso.{tag}.band_sigma",bp);   put(f"iso.{tag}.band_lam",bl)
    put(f"iso.{tag}.tau",bp*_col_s2cl2); put(f"iso.{tag}.short",TAU_REQ/(bp*_col_s2cl2))
    put(f"iso.{tag}.short_global",TAU_REQ/(gp*_col_s2cl2))
    put(f"iso.{tag}.ratio_global_band",gp/bp)
    pv,pl=_PUB[tag]
    put(f"iso.{tag}.val_pct",100*abs(gp-pv)/pv); put(f"iso.{tag}.val_nm",abs(gl-pl))

# five-point compatibility relation
_r=np.array([REG["e8b.L2pct"],REG["e11a.L2pct"],REG["e11b.L2pct"],
             REG["e8a.L2pct"],REG["e8c.L2pct"]])
_c=np.array([REG["corr.8B"],REG["corr.11A"],REG["corr.11B"],
             REG["corr.8A"],REG["corr.8C"]])
put("five.pearson", np.corrcoef(_r,_c)[0,1])
_rr=np.argsort(np.argsort(_r)).astype(float); _cc=np.argsort(np.argsort(_c)).astype(float)
put("five.spearman", np.corrcoef(_rr,_cc)[0,1])
put("iso.improve_sscl2", REG["scr.S2Cl2.short"]/REG["iso.SSCl2.short"])
put("iso.worsen_clsscl", REG["iso.ClSSCl.short"]/REG["scr.S2Cl2.short"])


# ============================================================================
# EXPERIMENTS 1 TO 6
# Experiment 1's reference is the pre-UUV laptop baseline, identified by
# reproducing all three published difference statistics exactly.
# ============================================================================
PRE = L("Nominal_Bkzz_SO2-25m0328-laptop-baseline.vul")
A05 = L("Experiment3_A0p5_UUV-half-25m0328.vul")
A2  = L("Experiment3_A2_UUV-double-25m0328.vul")
B4  = {b: L(f"Experiment4B_B{i}_{r}-25m0328.vul")
       for i, (b, r) in enumerate(
           [("B1","300-320"),("B2","320-340"),("B3","340-360"),
            ("B4","360-380"),("B5","380-420"),("B6","420-698")], start=1)}

def _diffstats(a, b):
    d = np.abs(np.asarray(a) - np.asarray(b))
    ref = np.abs(np.asarray(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(ref > 0, d / ref, 0.0)
    return float(d.max()), float(d.mean()), float(np.nanmean(rel))

# ---- Experiment 1: the A=0 case against the pre-UUV model -------------------
for key in ("ymix", "tau", "aflux"):
    if key in A0["variable"] and key in PRE["variable"]:
        mx, mn, mr = _diffstats(A0["variable"][key], PRE["variable"][key])
        put(f"e1.{key}.maxabs", mx); put(f"e1.{key}.meanabs", mn)
        put(f"e1.{key}.meanrel", mr)
_ja = np.array([A0["variable"]["J_sp"][k] for k in sorted(A0["variable"]["J_sp"])])
_jb = np.array([PRE["variable"]["J_sp"][k] for k in sorted(PRE["variable"]["J_sp"])])
mx, mn, mr = _diffstats(_ja, _jb)
put("e1.J_sp.maxabs", mx); put("e1.J_sp.meanabs", mn); put("e1.J_sp.meanrel", mr)

# ---- Experiments 2 and 3: species response against the A=0 control ----------
# Experiments 2 and 3 express relative change against the LARGER of the two
# values, not against the reference. That convention is nowhere stated in those
# documents; it was recovered by reproducing their published numbers, and it is
# the only definition that matches every row of the Experiment 2 table.
_sp = list(A0["variable"]["species"])
_m0 = np.asarray(A0["variable"]["ymix"])
SPECIES_1_6 = ("SO2", "SO3", "H2SO4", "H2SO4_l", "SO", "S2O2", "ClS2", "S2Cl2",
               "S2O", "S3", "SCl2", "SCl", "S4")
for tag, case in (("e2", A1), ("e3a", A05), ("e3b", A2)):
    _mc = np.asarray(case["variable"]["ymix"])
    for name in SPECIES_1_6:
        if name not in _sp:
            continue
        i = _sp.index(name)
        d = np.abs(_mc[:, i] - _m0[:, i])
        put(f"{tag}.{name}.maxabs", d.max())
        den = np.maximum(_m0[:, i], _mc[:, i])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(den > 0, d / den, 0.0)
        j = int(np.argmax(rel))
        put(f"{tag}.{name}.maxrel", rel[j])
        put(f"{tag}.{name}.level", j)
        put(f"{tag}.{name}.ref_at_level", _m0[j, i])
        put(f"{tag}.{name}.case_at_level", _mc[j, i])
        put(f"{tag}.{name}.signed_at_level", _mc[j, i] - _m0[j, i])

# ---- amplitude scaling, the four-point series ------------------------------
_f0 = np.asarray(A0["variable"]["aflux"])
for tag, case, amp in (("a05", A05, 0.5), ("a1", A1, 1.0), ("a2", A2, 2.0)):
    d = np.asarray(case["variable"]["aflux"]) - _f0
    put(f"amp.{tag}.L2", np.sqrt(np.sum(d ** 2 * quad[None, :])))
    put(f"amp.{tag}.maxabs", np.abs(d).max())
for a, b in (("a05", "a1"), ("a1", "a2"), ("a05", "a2")):
    put(f"amp.ratio_{b}_over_{a}", REG[f"amp.{b}.L2"] / REG[f"amp.{a}.L2"])

# ---- Experiment 4B: the six band runs --------------------------------------
_BANDS = {"B1": (300, 320), "B2": (320, 340), "B3": (340, 360),
          "B4": (360, 380), "B5": (380, 420), "B6": (420, 698)}
for b, case in B4.items():
    d = np.asarray(case["variable"]["aflux"]) - _f0
    put(f"e4b.{b}.L2", np.sqrt(np.sum(d ** 2 * quad[None, :])))
    put(f"e4b.{b}.maxabs", np.abs(d).max())
    ij = np.unravel_index(np.argmax(np.abs(d)), d.shape)
    put(f"e4b.{b}.peak_z", z[ij[0]]); put(f"e4b.{b}.peak_lam", bins[ij[1]])
    u = np.loadtxt(f"atm/UV_absorber_Experiment4B_{b}_25m0328.txt", skiprows=1)
    nz = u[u[:, 1] > 0]
    put(f"e4b.{b}.support_lo", nz[0, 0]); put(f"e4b.{b}.support_hi", nz[-1, 0])
    put(f"e4b.{b}.npoints", len(nz)); put(f"e4b.{b}.peakQ", nz[:, 1].max())
_tot = sum(REG[f"e4b.{b}.L2"] for b in B4)
for b in B4:
    put(f"e4b.{b}.share", 100 * REG[f"e4b.{b}.L2"] / _tot)

# ---------------------------------------------------------------------------
# Solver convergence diagnostics. The Experiment 1, 2 and 3 documents quote the
# run-termination state of each integration (longdy, longdydt, aflux_change and
# the step statistics) as evidence that the runs converged. These are stored
# scalars, not derived quantities: they are read straight out of each .vul.
# ---------------------------------------------------------------------------
LAP = L("Nominal_Bkzz_SO2-25m0328-laptop-baseline.vul")
NBL = L("Nominal_Bkzz_SO2-25m0328-baseline.vul")
NRF = L("Nominal_Bkzz_SO2-25m0328-reference.vul")
RUNS_ALL = {
    "e1.ref": LAP, "nom.base": NBL, "nom.ref": NRF,
    "e1.A0": A0, "e2.A1": A1, "e3a.A05": A05, "e3b.A2": A2,
    "e7a": E7A, "e7b": E7B, "e8.base": B8, "e8a": E8A, "e8b": E8B, "e8c": E8C,
    "e9": E9,
}
for b, case in B4.items():
    RUNS_ALL[f"e4b.{b}"] = case
for tag, sol in RUNS_ALL.items():
    for grp in ("variable", "parameter"):
        for k, v in sol[grp].items():
            if isinstance(v, (bool, str, dict)) or v is None:
                continue
            if np.isscalar(v) and np.isreal(v) and np.isfinite(float(v)):
                put(f"{tag}.{grp[0]}.{k}", float(v))

# ---------------------------------------------------------------------------
# Whole-array comparisons for Experiments 2 and 3.
#
# Two conventions matter here and both are taken from the documents themselves.
# (1) The mean is over EVERY element of the stored array, unweighted: no
#     quadrature weight and no altitude weight. It is a diagnostic of the saved
#     field, not a physical column integral.
# (2) The relative difference divides by max(|ref|, |case|), which the
#     Experiment 2 actinic-flux section prints explicitly as
#     max|dF / max(|F0|,|F1|)|. The abundance tables in Experiments 2 and 3 use
#     the same denominator, which is how that convention was recovered.
# ---------------------------------------------------------------------------
def array_pair(tag, field, ref, case):
    a = np.asarray(ref["variable"][field]); b = np.asarray(case["variable"][field])
    d = b - a
    ad = np.abs(d)
    put(f"{tag}.{field}.maxabs", ad.max())
    put(f"{tag}.{field}.meanabs", ad.mean())
    den = np.maximum(np.abs(a), np.abs(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(den > 0, ad / den, 0.0)
    put(f"{tag}.{field}.maxrel", rel.max())
    put(f"{tag}.{field}.meanrel", rel.mean())
    put(f"{tag}.{field}.medianrel", np.median(rel))
    ij = np.unravel_index(np.argmax(ad), ad.shape)
    put(f"{tag}.{field}.arg_i", ij[0]); put(f"{tag}.{field}.arg_j", ij[1])
    put(f"{tag}.{field}.ref_at_arg", a[ij]); put(f"{tag}.{field}.case_at_arg", b[ij])
    put(f"{tag}.{field}.delta_at_arg", d[ij])

def jsp_pair(tag, ref, case):
    Ja, Jb = ref["variable"]["J_sp"], case["variable"]["J_sp"]
    keys = sorted(set(Ja) & set(Jb), key=str)
    put(f"{tag}.J.nbranch", len(keys))
    put(f"{tag}.J.nvalues", sum(np.asarray(Ja[k]).size for k in keys))
    allmax = 0.0; tot = 0.0; n = 0
    for k in keys:
        a = np.asarray(Ja[k]); b = np.asarray(Jb[k]); ad = np.abs(b - a)
        allmax = max(allmax, ad.max()); tot += ad.sum(); n += ad.size
        den = np.maximum(np.abs(a), np.abs(b))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(den > 0, ad / den, 0.0)
        nm = f"{k[0]}_{k[1]}"
        put(f"{tag}.J.{nm}.maxabs", ad.max())
        put(f"{tag}.J.{nm}.maxrel", rel.max())
        put(f"{tag}.J.{nm}.meanrel", rel.mean())
    put(f"{tag}.J.maxabs", allmax); put(f"{tag}.J.meanabs", tot / n)
    _rs = []
    for k in keys:
        a = np.asarray(Ja[k]); b = np.asarray(Jb[k])
        den = np.maximum(np.abs(a), np.abs(b))
        with np.errstate(divide="ignore", invalid="ignore"):
            _rs.append(np.where(den > 0, np.abs(b - a) / den, 0.0))
    _rs = np.concatenate(_rs)
    put(f"{tag}.J.meanrel", _rs.mean()); put(f"{tag}.J.maxrel", _rs.max())

for tag, case in (("e2", A1), ("e3a", A05), ("e3b", A2)):
    for field in ("tau", "aflux", "sflux"):
        array_pair(tag, field, A0, case)
    jsp_pair(tag, A0, case)

# ---------------------------------------------------------------------------
# Experiment 1 compares its control run against the laptop baseline. That file is
# the reference identified by reproducing all three published diff statistics.
# ---------------------------------------------------------------------------
for field in ("tau", "aflux", "sflux"):
    array_pair("e1", field, LAP, A0)
jsp_pair("e1", LAP, A0)
_spL = list(LAP["variable"]["species"]); _mL = np.asarray(LAP["variable"]["ymix"])
_m0b = np.asarray(A0["variable"]["ymix"])
_dm = np.abs(_m0b - _mL)
put("e1.ymix.maxabs", _dm.max()); put("e1.ymix.meanabs", _dm.mean())
_den = np.maximum(np.abs(_mL), np.abs(_m0b))
with np.errstate(divide="ignore", invalid="ignore"):
    _rel = np.where(_den > 0, _dm / _den, 0.0)
put("e1.ymix.maxrel", _rel.max())
put("e1.ymix.meanrel", _rel.mean())
put("e1.ymix.medianrel", np.median(_rel))

# ---------------------------------------------------------------------------
# Conserved-element budgets. Experiment 3 tabulates the solver's own reported
# total atom loss per element, a stored dict, not a recomputed quantity.
# ---------------------------------------------------------------------------
for tag, sol in RUNS_ALL.items():
    for dn in ("atom_loss", "atom_sum", "atom_ini", "atom_conden"):
        if dn in sol["variable"]:
            for el, val in sol["variable"][dn].items():
                put(f"{tag}.{dn}.{el}", float(np.asarray(val)))

# ---------------------------------------------------------------------------
# Linearity of the response in A. If the perturbation were exactly linear the
# A=0.5 difference would be half the A=1 difference, so this ratio would be 1.
# The median is taken over the elements where the A=1 difference is non-zero,
# since the ratio is undefined where the denominator vanishes.
# ---------------------------------------------------------------------------
for field in ("tau", "aflux"):
    f0 = np.asarray(A0["variable"][field])
    fh = np.asarray(A05["variable"][field])
    f1 = np.asarray(A1["variable"][field])
    num = np.abs(fh - f0); den = 0.5 * np.abs(f1 - f0)
    m = den > 0
    put(f"e3.lin.{field}.median", np.median(num[m] / den[m]))
    put(f"e3.lin.{field}.mean", np.mean(num[m] / den[m]))

# The A=1 field sampled at the location where the A=0.5 difference peaks, which
# is how Experiment 3 tabulates all three cases side by side at one location.
for field in ("tau", "aflux"):
    f0 = np.asarray(A0["variable"][field]); fh = np.asarray(A05["variable"][field])
    f1 = np.asarray(A1["variable"][field]); f2 = np.asarray(A2["variable"][field])
    ij = np.unravel_index(np.argmax(np.abs(fh - f0)), f0.shape)
    put(f"e3a.{field}.arg_i", ij[0]); put(f"e3a.{field}.arg_j", ij[1])
    put(f"e3a.{field}.at_arg.A0", f0[ij]); put(f"e3a.{field}.at_arg.A05", fh[ij])
    put(f"e3a.{field}.at_arg.A1", f1[ij]); put(f"e3a.{field}.at_arg.A2", f2[ij])

# ---------------------------------------------------------------------------
# Linearity of the response in A, restricted to where the absorber acts.
#
# The absorber carries opacity only on 300 to 698 nm. Outside that range the
# A=0.5 and A=1 fields differ from the control only by solver round-off: the
# median |difference| there is 5e-7 of the field itself, against 1e-3 inside.
# A ratio of two round-off differences carries no information, and the out-of-band
# elements are 83790 of the 95190 in the optical-depth array, so any statistic
# taken over the whole array is dominated by noise. The ratio is therefore
# reported over the absorber's own wavelength range, where both differences are
# real signal. Under that definition the response is linear in A.
# ---------------------------------------------------------------------------
_inb = (bins >= 300.0) & (bins <= 698.0)
for field in ("tau", "aflux"):
    f0 = np.asarray(A0["variable"][field]); fh = np.asarray(A05["variable"][field])
    f1 = np.asarray(A1["variable"][field])
    num = np.abs(fh - f0); den = 0.5 * np.abs(f1 - f0)
    for nm, w in (("inband", _inb), ("outband", ~_inb), ("all", np.ones_like(_inb))):
        m = (den > 0) & w[None, :]
        r = num[m] / den[m]
        put(f"e3.lin.{field}.{nm}.median", np.median(r))
        put(f"e3.lin.{field}.{nm}.n", int(m.sum()))
        sc = np.abs(f0)[m]
        put(f"e3.lin.{field}.{nm}.median_relsize",
            np.median(np.where(sc > 0, num[m] / np.maximum(sc, 1e-300), 0.0)))

# ---------------------------------------------------------------------------
# Experiment 5: exact tendency decomposition.
#
# The generated chem_funs.py builds each species' chemical tendency as a sum of
# signed reaction terms. This block parses those terms out of the generated
# source, evaluates each at a saved state, and checks they sum to the matching
# column of the full chemdf() result. Summing in source order reproduces the
# model's own arithmetic exactly, so the reconstruction error is identically
# zero rather than the 1e-7 round-off an out-of-order summation leaves.
#
# Terms sharing a reaction id are added together before ranking: a reaction that
# appears twice in one equation contributes once to that species' tendency.
# ---------------------------------------------------------------------------
sys.path.insert(0, ".")
import re as _re
from collections import defaultdict as _dd
import chem_funs as _cf

_SRC = open("chem_funs.py").read()
_BODY = _re.search(r"def chemdf\(y, M, k\):.*?\n(?=def )", _SRC, _re.S).group(0)

def _terms(sp_index):
    line = _re.search(rf"^\s*dydt\[{sp_index}\]\s*=(.*)$", _BODY, _re.M).group(1)
    return [(1 if m.group(1) == "+" else -1, int(m.group(3)), m.group(2), m.group(4))
            for m in _re.finditer(r"([+-])1\*(v_(\d+))\(([^()]*)\)", line)]

def _contrib(sol, sp_name):
    var, atm = sol["variable"], sol["atm"]
    i = list(var["species"]).index(sp_name)
    yT = np.transpose(np.asarray(var["y"]))
    M = np.asarray(atm["M"]); k = var["k"]
    rows, ids = [], []
    for sign, rid, fname, args in _terms(i):
        vals = eval(f"[{args}]", {"k": k, "M": M, "y": yT, "np": np})
        rows.append(sign * np.asarray(getattr(_cf, fname)(*vals), dtype=float))
        ids.append(rid)
    return np.array(rows), ids, i

def _chemdf(sol):
    return np.asarray(_cf.chemdf(np.asarray(sol["variable"]["y"]),
                                 np.asarray(sol["atm"]["M"]), sol["variable"]["k"]))

E5_SPECIES = ("SO2", "SO", "SO3", "H2SO4", "H2SO4_l")
for sp in E5_SPECIES:
    agg = {}
    for nm, sol in (("A0", A0), ("A1", A1)):
        rows, ids, i = _contrib(sol, sp)
        recon = rows.sum(axis=0)
        put(f"e5.{sp}.{nm}.extraction_error", np.abs(recon - _chemdf(sol)[:, i]).max())
        put(f"e5.{sp}.nterms", len(ids)); put(f"e5.{sp}.nunique", len(set(ids)))
        a = _dd(lambda: np.zeros(rows.shape[1]))
        for rid, row in zip(ids, rows):
            a[rid] = a[rid] + row
        agg[nm] = a
    for rid in agg["A0"]:
        d = np.abs(agg["A1"][rid] - agg["A0"][rid])
        j = int(np.argmax(d))
        put(f"e5.{sp}.R{rid}.maxdelta", d[j])
        put(f"e5.{sp}.R{rid}.level", j)
        put(f"e5.{sp}.R{rid}.bg_A0", agg["A0"][rid][j])
        put(f"e5.{sp}.R{rid}.bg_A1", agg["A1"][rid][j])
        put(f"e5.{sp}.R{rid}.signed_delta", agg["A1"][rid][j] - agg["A0"][rid][j])
        if agg["A0"][rid][j] != 0:
            put(f"e5.{sp}.R{rid}.frac_pct", 100 * d[j] / abs(agg["A0"][rid][j]))
            put(f"e5.{sp}.R{rid}.frac", d[j] / abs(agg["A0"][rid][j]))

# ---------------------------------------------------------------------------
# Experiment 6: where in altitude the response sits.
#
# The vertical weight is the quadrature-weighted absolute actinic-flux change
# summed over ALL wavelength bins, not only the band the run perturbed. That is
# the definition that reproduces the published per-band mean altitudes, and it is
# the physically sensible one: a band-limited opacity change alters the radiation
# field outside its own window too, and that response belongs in the total.
# ---------------------------------------------------------------------------
def vweight(case, ref=A0):
    d = np.abs(np.asarray(case["variable"]["aflux"]) - np.asarray(ref["variable"]["aflux"]))
    return np.maximum((d * quad[None, :]).sum(axis=1), 0.0)

def wquantiles(w, zz=z):
    tot = w.sum(); c = np.cumsum(w) / tot
    return [float(zz[np.searchsorted(c, p)]) for p in (0.25, 0.50, 0.75)]

def narrowest(w, zz, frac):
    """Narrowest contiguous altitude interval holding at least frac of the signal."""
    tot = w.sum(); n = len(w); best = None
    for i in range(n):
        acc = 0.0
        for j in range(i, n):
            acc += w[j]
            if acc >= frac * tot:
                width = zz[j] - zz[i]
                if best is None or width < best[0]:
                    best = (width, zz[i], zz[j], 100.0 * acc / tot)
                break
    return best

for b, case in B4.items():
    w = vweight(case)
    put(f"e6.{b}.meanz", np.sum(z * w) / w.sum())
    q = wquantiles(w)
    for nm, v in zip(("q25", "q50", "q75"), q):
        put(f"e6.{b}.{nm}", v)

_wtot = vweight(A1)
put("e6.A1.meanz", np.sum(z * _wtot) / _wtot.sum())
for nm, v in zip(("q25", "q50", "q75"), wquantiles(_wtot)):
    put(f"e6.A1.{nm}", v)
for frac in (0.50, 0.68, 0.80, 0.90, 0.95, 0.99):
    r = narrowest(_wtot, z, frac)
    if r:
        put(f"e6.interval.{int(frac*100)}.width", r[0])
        put(f"e6.interval.{int(frac*100)}.lo", r[1])
        put(f"e6.interval.{int(frac*100)}.hi", r[2])
        put(f"e6.interval.{int(frac*100)}.contained", r[3])

# Containment of each reaction's response inside the response band.
#
# The quantity is the reaction's signed contribution to a species tendency, and
# the containment is the share of the total absolute change that falls in
# 69 to 93 km. It is species-independent: a reaction enters two equations with
# opposite signs but the same magnitude, so reading R781 out of the SO2 equation
# and out of the SO equation gives the same 99.1283 per cent. That is why this is
# reported per reaction rather than per species.
_aggs = {}
for sp in E5_SPECIES:
    for nm, sol in (("A0", A0), ("A1", A1)):
        rows, ids, _ = _contrib(sol, sp)
        a = _dd(lambda: np.zeros(rows.shape[1]))
        for rid, row in zip(ids, rows):
            a[rid] = a[rid] + row
        _aggs[(sp, nm)] = a
for sp in E5_SPECIES:
    a0, a1 = _aggs[(sp, "A0")], _aggs[(sp, "A1")]
    for rid in a0:
        dd = np.abs(a1[rid] - a0[rid])
        if dd.sum() <= 0:
            continue
        # Experiment 6 partitions the column into three regions and reports the
        # share of each reaction's response in each: below the cloud base, the
        # intermediate layer, and the response band itself.
        # The three regions must partition the column. 59 km is a grid level, so it
        # must belong to exactly one of them: it is counted in the 59 to 68 km
        # layer, not in "below 59". Including it in both makes the shares sum to
        # more than 100 per cent.
        for _bn, _bm in (("below59", z < 59.0), ("b59_68", (z >= 59.0) & (z <= 68.0)),
                         ("band69_93", band)):
            put(f"e6.R{rid}.pct_{_bn}", 100.0 * dd[_bm].sum() / dd.sum())
        put(f"e6.R{rid}.pct_in_band", 100.0 * dd[band].sum() / dd.sum())
        put(f"e6.R{rid}.peak_z", z[int(np.argmax(dd))])
        put(f"e6.R{rid}.peak_val", dd.max())

# The photolysis rate coefficients themselves, for the Experiment 6 J tables.
_Ja, _Jb = A0["variable"]["J_sp"], A1["variable"]["J_sp"]
for kk in sorted(set(_Ja) & set(_Jb), key=str):
    dd = np.abs(np.asarray(_Jb[kk]) - np.asarray(_Ja[kk]))
    if dd.sum() <= 0:
        continue
    nm = f"{kk[0]}_{kk[1]}"
    put(f"e6.J.{nm}.pct_in_band", 100.0 * dd[band].sum() / dd.sum())
    put(f"e6.J.{nm}.peak_z", z[int(np.argmax(dd))])

# ---------------------------------------------------------------------------
# Experiment 4, Stage 4A: the exact spectral kernel.
#
# VULCAN integrates each photolysis branch with a trapezoidal rule that uses two
# bin widths, 0.1 nm below the transition bin and 2 nm above it, and subtracts
# half of each segment's two endpoint values. Reproducing that rule exactly gives
# back the stored J_sp with an error of identically zero, which is what makes the
# band decomposition below the model's own arithmetic rather than a separate
# numerical scheme.
#
# Two conventions to keep straight:
#  - A band integral applies the same endpoint half-weighting at the band edges.
#    Omitting it inflates every band by a few per cent.
#  - The fractions divide by the STORED full-spectrum change, not by the sum of
#    the band contributions, so a branch's fractions need not total exactly 100.
# ---------------------------------------------------------------------------
IDX = int(np.searchsorted(bins, 240.0))          # 1440; bins[1440] is 240.0 nm
DB1 = float(A1["variable"]["dbin1"]); DB2 = float(A1["variable"]["dbin2"])
WQ = np.where(bins < 240.0, DB1, DB2)
E4A_BANDS = [(300, 320), (320, 340), (340, 360), (360, 380), (380, 420), (420, 620)]
LEVEL_4A = 41

def vulcan_J(flux, sigma):
    """VULCAN's own photolysis quadrature, reproduced exactly."""
    r = np.sum(flux[:, :IDX] * sigma[None, :IDX] * DB1, axis=1)
    r -= 0.5 * (flux[:, 0] * sigma[0] + flux[:, IDX - 1] * sigma[IDX - 1]) * DB1
    r += np.sum(flux[:, IDX:] * sigma[None, IDX:] * DB2, axis=1)
    r -= 0.5 * (flux[:, IDX] * sigma[IDX] + flux[:, -1] * sigma[-1]) * DB2
    return r

def band_integral(row, sigma, lo, hi):
    m = np.where((bins >= lo) & (bins <= hi))[0]
    r = np.sum(row[m] * sigma[m] * WQ[m])
    r -= 0.5 * (row[m[0]] * sigma[m[0]] + row[m[-1]] * sigma[m[-1]]) * WQ[m[0]]
    return r

_cJ = A1["variable"]["cross_J"]
_J0, _J1 = A0["variable"]["J_sp"], A1["variable"]["J_sp"]
_dF = np.asarray(A1["variable"]["aflux"]) - np.asarray(A0["variable"]["aflux"])

for kk in sorted(_cJ, key=str):
    if kk not in _J1:
        continue
    sig = np.asarray(_cJ[kk]); nm = f"{kk[0]}_{kk[1]}"
    # the reconstruction check, on both states
    for tg, sol in (("A0", A0), ("A1", A1)):
        err = np.abs(vulcan_J(np.asarray(sol["variable"]["aflux"]), sig)
                     - np.asarray(sol["variable"]["J_sp"][kk]))
        put(f"e4a.{nm}.{tg}.recon_maxerr", err.max())
        put(f"e4a.{nm}.{tg}.recon_meanerr", err.mean())
    den = _J1[kk][LEVEL_4A] - _J0[kk][LEVEL_4A]
    raw = np.array([band_integral(_dF[LEVEL_4A], sig, lo, hi) for lo, hi in E4A_BANDS])
    for (lo, hi), rv in zip(E4A_BANDS, raw):
        put(f"e4a.{nm}.{lo}_{hi}.dJ", rv)
        if den != 0:
            put(f"e4a.{nm}.{lo}_{hi}.pct", 100.0 * rv / den)

# Combined L2 of all explicit branch-rate changes, per Stage 4A band.
for lo, hi in E4A_BANDS:
    tot = 0.0
    for kk in _cJ:
        if kk not in _J1 or kk[1] == 0:
            continue
        sig = np.asarray(_cJ[kk])
        col = np.array([band_integral(_dF[i], sig, lo, hi) for i in range(_dF.shape[0])])
        tot += np.sum(col ** 2)
    put(f"e4a.L2.{lo}_{hi}", np.sqrt(tot))

# ---------------------------------------------------------------------------
# Experiment 4, Stage 4B: isolated-window response ratios.
#
# Each of the six windows is its own steady-state run, so the ratio compares that
# run's stored photolysis rate against the control, over the full-spectrum change:
#     R_B = (J_B - J_A0) / (J_A1 - J_A0)   at level 41.
# These are cross-run comparisons, not a decomposition, so they need not sum to
# 100 per cent, and Stage 4A's fractions and these ratios are independent evidence.
# ---------------------------------------------------------------------------
for b, case in B4.items():
    Jb = case["variable"]["J_sp"]
    for kk in sorted(set(Jb) & set(_J1), key=str):
        den = _J1[kk][LEVEL_4A] - _J0[kk][LEVEL_4A]
        if den == 0:
            continue
        num = Jb[kk][LEVEL_4A] - _J0[kk][LEVEL_4A]
        nm = f"{kk[0]}_{kk[1]}"
        put(f"e4b.{b}.{nm}.ratio_pct", 100.0 * num / den)
        put(f"e4b.{b}.{nm}.dJ", num)

# Steady-state sulfur response of each isolated window against the control.
_sp0 = list(A0["variable"]["species"]); _ym0 = np.asarray(A0["variable"]["ymix"])
for b, case in B4.items():
    _ymb = np.asarray(case["variable"]["ymix"])
    for name in SPECIES_1_6:
        if name not in _sp0:
            continue
        i = _sp0.index(name)
        d = np.abs(_ymb[:, i] - _ym0[:, i])
        put(f"e4b.{b}.{name}.maxabs", d.max())
        put(f"e4b.{b}.{name}.level", int(np.argmax(d)))
        den = np.maximum(_ym0[:, i], _ymb[:, i])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(den > 0, d / den, 0.0)
        put(f"e4b.{b}.{name}.maxrel", rel.max())

# ---------------------------------------------------------------------------
# Experiment 3, four-point amplitude sensitivity (A = 0, 0.5, 1, 2).
#
# Linearity is tested on the L2 norm of the response field, not element by
# element: the norm is carried by the elements where the absorber actually acts,
# so it measures the physical response rather than round-off. The residual is the
# norm of the departure from exact proportionality, ||dF(A) - A dF(1)|| / ||A dF(1)||,
# which is why it is larger than |1 - ratio|.
# ---------------------------------------------------------------------------
CLOUD_BASE_LEVEL = int(np.argmin(np.abs(z - 59.0)))   # level 30, z = 59.0 km
put("e3fp.cloud_base_level", CLOUD_BASE_LEVEL)
put("e3fp.cloud_base_z", z[CLOUD_BASE_LEVEL])
_AF = {a: np.asarray(s_["variable"]["aflux"]) - np.asarray(A0["variable"]["aflux"])
       for a, s_ in ((0.5, A05), (1.0, A1), (2.0, A2))}
_TAUF = {a: np.asarray(s_["variable"]["tau"]) - np.asarray(A0["variable"]["tau"])
         for a, s_ in ((0.5, A05), (1.0, A1), (2.0, A2))}
def _l2(x, w=None):
    return np.sqrt(np.sum(x ** 2 if w is None else x ** 2 * w))
for fld, DD in (("aflux", _AF), ("tau", _TAUF)):
    put(f"e3fp.{fld}.norm.A1", _l2(DD[1.0]))
    for a in (0.5, 2.0):
        tagA = str(a).replace(".", "p")
        put(f"e3fp.{fld}.norm.A{tagA}", _l2(DD[a]))
        put(f"e3fp.{fld}.ratio.A{tagA}", _l2(DD[a]) / (a * _l2(DD[1.0])))
        put(f"e3fp.{fld}.resid.A{tagA}", _l2(DD[a] - a * DD[1.0]) / _l2(a * DD[1.0]))
        put(f"e3fp.{fld}.resid_pct.A{tagA}",
            100.0 * _l2(DD[a] - a * DD[1.0]) / _l2(a * DD[1.0]))

# Four-point chemical amplitude response, per species, against the control.
for tag, case in (("A0p5", A05), ("A1", A1), ("A2", A2)):
    _ymc = np.asarray(case["variable"]["ymix"])
    for name in SPECIES_1_6:
        if name not in _sp:
            continue
        i = _sp.index(name)
        d = np.abs(_ymc[:, i] - _m0[:, i])
        j = int(np.argmax(d))
        put(f"e3fp.{tag}.{name}.maxabs", d.max())
        put(f"e3fp.{tag}.{name}.level", j)
        den = np.maximum(_m0[:, i], _ymc[:, i])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(den > 0, d / den, 0.0)
        k = int(np.argmax(rel))
        put(f"e3fp.{tag}.{name}.maxrel", rel[k])
        put(f"e3fp.{tag}.{name}.maxrel_pct", 100.0 * rel[k])
        put(f"e3fp.{tag}.{name}.maxrel_level", k)
        put(f"e3fp.{tag}.{name}.ref_at", _m0[k, i])
        put(f"e3fp.{tag}.{name}.case_at", _ymc[k, i])
        put(f"e3fp.{tag}.{name}.signed_at", _ymc[k, i] - _m0[k, i])
# amplitude ratios of the per-species response
for name in SPECIES_1_6:
    if name not in _sp:
        continue
    b = REG.get(f"e3fp.A1.{name}.maxabs")
    if not b:
        continue
    for tag, a in (("A0p5", 0.5), ("A2", 2.0)):
        v = REG.get(f"e3fp.{tag}.{name}.maxabs")
        if v:
            put(f"e3fp.{tag}.{name}.amp_ratio", v / (a * b))

# ---------------------------------------------------------------------------
# Bounded relative change of the stored fields, per amplitude.
#
# "Bounded" means the denominator is max(|ref|, |case|), so the ratio can never
# exceed 1. The mean is taken over the elements where that denominator is
# non-zero: an element where both states are exactly zero has no relative change
# to average in, and including it as a zero would dilute the mean by an arbitrary
# amount depending on how much of the grid is dark. This definition reproduces
# the optical-depth table exactly.
# ---------------------------------------------------------------------------
for tag, case in (("A0p5", A05), ("A1", A1), ("A2", A2)):
    for fld in ("tau", "aflux", "sflux"):
        a = np.asarray(A0["variable"][fld]); b = np.asarray(case["variable"][fld])
        ad = np.abs(b - a); den = np.maximum(np.abs(a), np.abs(b))
        m = den > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(m, ad / den, 0.0)
        put(f"e3fp.bnd.{tag}.{fld}.max", rel[m].max())
        put(f"e3fp.bnd.{tag}.{fld}.max_pct", 100.0 * rel[m].max())
        put(f"e3fp.bnd.{tag}.{fld}.mean", rel[m].mean())
        put(f"e3fp.bnd.{tag}.{fld}.mean_pct", 100.0 * rel[m].mean())
        put(f"e3fp.bnd.{tag}.{fld}.n", int(m.sum()))
        _mr = m[CLOUD_BASE_LEVEL] if rel.ndim == 2 else None
        if rel.ndim == 2 and m[CLOUD_BASE_LEVEL].any():
            _r = rel[CLOUD_BASE_LEVEL][m[CLOUD_BASE_LEVEL]]
            put(f"e3fp.bnd.{tag}.{fld}.cb_max_pct", 100.0 * _r.max())
            put(f"e3fp.bnd.{tag}.{fld}.cb_mean_pct", 100.0 * _r.mean())

# Same bounded measure per photolysis branch, which is how the branch table ranks.
for tag, case in (("A0p5", A05), ("A1", A1), ("A2", A2)):
    Jc = case["variable"]["J_sp"]
    for kk in sorted(set(_J0) & set(Jc), key=str):
        a = np.asarray(_J0[kk]); b = np.asarray(Jc[kk])
        ad = np.abs(b - a); den = np.maximum(np.abs(a), np.abs(b))
        m = den > 0
        if not m.any():
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(m, ad / den, 0.0)
        nm = f"{kk[0]}_{kk[1]}"
        put(f"e3fp.bnd.{tag}.J.{nm}.max_pct", 100.0 * rel[m].max())
        put(f"e3fp.bnd.{tag}.J.{nm}.mean_pct", 100.0 * rel[m].mean())
        put(f"e3fp.bnd.{tag}.J.{nm}.argmax_level", int(np.argmax(rel)))
        # The four-point branch tables report the response at the cloud base,
        # level 30 = 59.0 km, which is the same boundary Experiment 6 partitions
        # on. They are captioned as maxima over altitude, but every published
        # value is the level-30 value; the true maximum sits at 53 to 55 km,
        # inside the lower cloud, and is roughly half again as large.
        put(f"e3fp.bnd.{tag}.J.{nm}.at_cloudbase_pct", 100.0 * rel[CLOUD_BASE_LEVEL])

# L2 amplitude test on the sulfur reservoirs themselves. Unlike the radiative
# response, these ratios sit far from 1, so the chemical response to amplitude is
# not proportional even though the flux response is.
for name in SPECIES_1_6:
    if name not in _sp:
        continue
    i = _sp.index(name)
    d1 = np.asarray(A1["variable"]["ymix"])[:, i] - _m0[:, i]
    n1 = np.sqrt(np.sum(d1 ** 2))
    if n1 == 0:
        continue
    for tag, a, case in (("A0p5", 0.5, A05), ("A2", 2.0, A2)):
        da = np.asarray(case["variable"]["ymix"])[:, i] - _m0[:, i]
        put(f"e3fp.L2.{tag}.{name}.ratio", np.sqrt(np.sum(da ** 2)) / (a * n1))
        put(f"e3fp.L2.{tag}.{name}.norm", np.sqrt(np.sum(da ** 2)))
    put(f"e3fp.L2.A1.{name}.norm", n1)

# ---------------------------------------------------------------------------
# Experiment 7: reaction-tendency response to vertical displacement.
#
# Units matter here and the published table got them mixed. The tendency is a
# volumetric rate, cm^-3 s^-1. Integrating it over altitude gives a column rate,
# cm^-2 s^-1, only if the layer thickness is in cm. The published integrated
# column used the thickness in km, which leaves the value 1e5 too small and the
# units mixed. Both forms are registered: ".col" is the correct column integral
# in cm^-2 s^-1, ".col_km" is the quantity as printed. The mantissas are
# identical, so only the exponent distinguishes them.
# ---------------------------------------------------------------------------
E7_LEVEL_71 = int(np.argmin(np.abs(z - 71.0)))
E7_LEVEL_69 = int(np.argmin(np.abs(z - 69.0)))
put("e7.level71", E7_LEVEL_71); put("e7.level69", E7_LEVEL_69)

def _agg_by_reaction(sol, sp):
    rows, ids, _ = _contrib(sol, sp)
    a = _dd(lambda: np.zeros(rows.shape[1]))
    for rid, row in zip(ids, rows):
        a[rid] = a[rid] + row
    return a

for sp in E5_SPECIES + ("SO",):
    if sp not in _sp:
        continue
    base = _agg_by_reaction(A1, sp)
    for tag, case in (("7a", E7A), ("7b", E7B)):
        cur = _agg_by_reaction(case, sp)
        for rid in base:
            d = cur[rid] - base[rid]
            if not np.any(d):
                continue
            k = f"e7.{tag}.{sp}.R{rid}"
            put(f"{k}.col", np.sum(d[band] * dzcm[band]))
            put(f"{k}.col_km", np.sum(d[band] * dzkm[band]))
            put(f"{k}.at71", d[E7_LEVEL_71])
            put(f"{k}.at69", d[E7_LEVEL_69])
            put(f"{k}.peakabs", np.abs(d).max())

# Signed global integral of the actinic-flux change, quadrature weighted.
for tag, case in (("7a", E7A), ("7b", E7B)):
    d = np.asarray(case["variable"]["aflux"]) - np.asarray(A1["variable"]["aflux"])
    put(f"e7.{tag}.signed_global", np.sum(d * quad[None, :]))
    put(f"e7.{tag}.signed_global_plain", d.sum())

# Experiment 7 radiative norms, absorber columns and realized centroids.
#
# The displacement is NOT symmetric. Shifting the profile up truncates its tail at
# the top of the grid, so renormalization leaves the realized centroid 17.844 km
# above nominal, while shifting down has room and lands exactly 18 km below. The
# published table reports both as +/- 17.844 km.
_F0 = np.asarray(A1["variable"]["aflux"])
put("e7.nominal_norm", np.sqrt(np.sum(_F0 ** 2 * quad[None, :])))
for tag, case in (("7a", E7A), ("7b", E7B)):
    d = np.asarray(case["variable"]["aflux"]) - _F0
    put(f"e7.{tag}.L2_quad", np.sqrt(np.sum(d ** 2 * quad[None, :])))
    put(f"e7.{tag}.L2_plain", np.sqrt(np.sum(d ** 2)))
    put(f"e7.{tag}.L2_share_pct", 100 * np.sqrt(np.sum(d ** 2 * quad[None, :])) / REG["e7.nominal_norm"])

import glob as _glob
_PROFS = {"nominal": "atm/mode1+2.txt",
          "7a": "atm/mode1+2_Experiment7A_shift18_25m0328.txt",
          "7b": "atm/mode1+2_Experiment7B_shiftDown18_25m0328.txt"}
_cols = {}
for tag, pat in _PROFS.items():
    g = _glob.glob(pat)
    if not g:
        continue
    N = Nprof(g[0]); col = float(np.sum(N * dzcm))
    _cols[tag] = col
    put(f"e7.{tag}.column", col)
    put(f"e7.{tag}.centroid", np.sum(zu * N * dzcm) / col)
if "nominal" in _cols:
    for tag in ("7a", "7b"):
        if tag in _cols:
            put(f"e7.{tag}.centroid_shift",
                REG[f"e7.{tag}.centroid"] - REG["e7.nominal.centroid"])
            put(f"e7.{tag}.column_reldiff",
                abs(_cols[tag] - _cols["nominal"]) / _cols["nominal"])

# Share of the integrated model Q-weight in the 620 to 698 nm tail. The boundary
# convention matters at the second decimal: 620 nm is the last tabulated entry, so
# excluding it gives 0.642 per cent and including it gives 0.664 per cent.
_Qg = interpolate.interp1d(np.loadtxt("atm/UV_absorber.txt", skiprows=1)[:, 0],
                           np.loadtxt("atm/UV_absorber.txt", skiprows=1)[:, 1],
                           bounds_error=False, fill_value=0.0)(bins)
_all = (bins >= 300) & (bins <= 698)
_tot = np.sum(_Qg[_all] * quad[_all])
for nm, m in (("excl", (bins > 620) & (bins <= 698)), ("incl", (bins >= 620) & (bins <= 698))):
    put(f"e7.tail_share_pct.{nm}", 100 * np.sum(_Qg[m] * quad[m]) / _tot)

# ---------------------------------------------------------------------------
# Experiment 13: Candidate A against the published bulk-absorbance requirement.
#
# Spacek et al. (2026) derive a decadic absorbance of 1278 cm^-1 at 375 nm for the
# bulk liquid of a Venus cloud droplet. Jiang et al. (2024) give mass extinction
# coefficients for the two iron-sulfur minerals. Dividing one by the other gives
# the mass loading each mineral would need, which is then compared with the droplet
# mass implied by the stated solution density.
#
# Two conservatisms are worth recording. Jiang's coefficient is an EXTINCTION
# coefficient measured through a suspension, so it includes scattering and
# therefore overstates absorption, which makes every shortfall here a lower bound.
# And the requirement is decadic; reading Jiang's optical depth as napierian
# instead would raise the required loading by ln(10), not lower it.
# ---------------------------------------------------------------------------
A_REQ_DECADIC = 1278.0        # cm^-1 at 375 nm, Spacek et al. (2026)
RHO_SOLUTION = 1.8            # g cm^-3, the droplet solution density they state
LAMBDA_REQ = 375.0
put("e13.A_required", A_REQ_DECADIC)
put("e13.rho_solution", RHO_SOLUTION)

_E13 = os.path.join(os.path.dirname(os.path.abspath(__file__)))
_MIN = {"rhomboclase": "jiang_rhomboclase_extracted.dat",
        "acid_ferric_sulfate": "jiang_acid-ferric-sulfate_extracted.dat"}
for _nm, _fn in _MIN.items():
    _cands = [os.path.join("Experiment13", "data", _fn), os.path.join(_E13, _fn)]
    _path = next((c for c in _cands if os.path.exists(c)), None)
    if _path is None:
        continue
    d = np.loadtxt(_path)
    lam, eps = d[:, 0], d[:, 1]
    put(f"e13.{_nm}.peak_eps", eps.max())
    put(f"e13.{_nm}.peak_lam", lam[int(np.argmax(eps))])
    e375 = float(interpolate.interp1d(lam, eps)(LAMBDA_REQ))
    put(f"e13.{_nm}.eps_375", e375)
    # required loading in g cm^-3, then kg per litre (numerically the same number)
    load = A_REQ_DECADIC / e375
    put(f"e13.{_nm}.loading_kg_per_L", load)
    put(f"e13.{_nm}.pct_of_droplet", 100.0 * load / RHO_SOLUTION)
    # the most favourable case: at the mineral's own peak
    put(f"e13.{_nm}.loading_at_peak", A_REQ_DECADIC / eps.max())
    put(f"e13.{_nm}.pct_at_peak", 100.0 * (A_REQ_DECADIC / eps.max()) / RHO_SOLUTION)
    # if Jiang's optical depth were napierian rather than decadic
    put(f"e13.{_nm}.loading_napierian", load * np.log(10.0))

# Method calibration: reproducing the paper's own ferric chloride figure. They quote
# 1 absorbance unit per g/L, so the requirement implies 1.278 kg/L against their 1.3.
put("e13.fecl3.eps_per_gpl", 1.0)
put("e13.fecl3.loading_kg_per_L", A_REQ_DECADIC / (1.0 / 0.001))
put("e13.fecl3.published_loading", 1.3)

# ---------------------------------------------------------------------------
# Experiment 12: the three-axis admissibility screen.
#
# The screen's arithmetic is not reimplemented here. The screen writes its own
# results to Experiment12/screen_results_25m0328.json and this block folds that
# file into the registry. One implementation, so the two cannot drift apart, which
# is the failure mode that produced most of the defects this harness has found.
#
# If the JSON is absent, run:
#     python Experiment12/three_axis_screen_25m0328.py
# ---------------------------------------------------------------------------
_E12 = os.path.join("Experiment12", "screen_results_25m0328.json")
if os.path.exists(_E12):
    with open(_E12) as _fh:
        _sc = json.load(_fh)
    for _k in ("column_aloft", "n_screened", "n_tracked", "n_nosupport", "n_nocross",
               "n_cov480", "n_cov400", "n_gap", "best_shortfall", "best_full_shortfall",
               "second_full_shortfall", "full_separation"):
        if _k in _sc and _sc[_k] is not None:
            put(f"e12.{_k}", _sc[_k])
    for _l, _v in _sc.get("q_at_probe", {}).items():
        put(f"e12.q.{_l}", _v)
    for _l, _v in _sc.get("tau_required_aloft", {}).items():
        put(f"e12.tau_req_aloft.{_l}", _v)
    for _name, _d in _sc.get("species", {}).items():
        for _f in ("col58", "col69", "centroid", "frac_in_band_pct", "peak_lam",
                   "peak_sigma", "tau58", "tau69", "shortfall58", "corr",
                   "share_320_400_pct", "support_lo", "support_hi"):
            if _d.get(_f) is not None:
                put(f"e12.{_name}.{_f}", _d[_f])
        for _l, _v in _d.get("shortfall_at", {}).items():
            if _v is not None:
                put(f"e12.{_name}.short.{_l}", _v)
else:
    print(f"NOTE: {_E12} not found, so Experiment 12 is absent from the registry.")
    print("      Run: python Experiment12/three_axis_screen_25m0328.py")

REG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "registry.json")
json.dump(REG, open(REG_PATH, "w"), indent=1)
print(f"registry written to {REG_PATH}")
print(f"registry entries: {len(REG)}")
for k in sorted(REG):
    print(f"  {k:32s} {REG[k]:.6g}")
