#!/usr/bin/env python3
"""Cross-check numeric literals in the .tex documents against the recomputed registry.
A literal matches a registry value w if the literal equals w rounded to the
literal's own number of significant digits, within half a unit in the last place."""
import re, json, sys, math, os
REG=json.load(open(os.path.join(os.path.dirname(__file__),"registry.json")))
vals=[(k,w) for k,w in REG.items() if math.isfinite(w) and w!=0]

def sig_of(t):
    t=t.lstrip('+-')
    if '.' not in t: return len(t.lstrip('0')) or 1
    a,b=t.split('.')
    a=a.lstrip('0')
    return len(a)+len(b) if a else max(len(b.lstrip('0')),1)

def match(v,s):
    av=abs(v)
    if av==0: return []
    out=[]
    for k,w in vals:
        aw=abs(w)
        e=math.floor(math.log10(aw))
        ulp=10**(e-s+1)
        # Never demand better agreement than double precision can express. A 16
        # digit literal implies a tolerance below one double ulp, and the literal
        # itself is reconstructed as mantissa*10**exp, which rounds. Floor the
        # tolerance at a few ulps of the value being compared.
        tol=max(0.5*ulp*1.001, 4*math.ulp(aw))
        if abs(av-aw) <= tol: out.append(k)
    return out

# The exponent may or may not be braced: both \times10^5 and \times10^{-5} occur
# across these documents, and an unbraced exponent previously fell through to the
# plain-decimal rule, which then compared the bare mantissa and reported a false
# mismatch. \cdot is also used in a few places.
SCI=re.compile(r'([+-]?\d+\.?\d*)\s*(?:\\times|\\cdot)\s*10\^\{?(-?\d+)\}?')
PCT=re.compile(r'([+-]?\d+\.\d{2,})\s*\\%')
DEC=re.compile(r'(?<![\d.])(\d+\.\d{3,})(?!\d)')

def scan(path):
    s=re.sub(r'(?m)^%.*','',open(path).read())
    # LaTeX lengths and font sizes are typography, not data
    s=re.sub(r'\{[\d.]+\\(?:text|line|column|page)width\}', ' ', s)
    s=re.sub(r'\\(?:setlength|hspace|vspace|tabcolsep|emergencystretch)[^\n]*', ' ', s); out=[]
    def eat(rx,conv,pct=False):
        nonlocal s
        for m in list(rx.finditer(s)):
            out.append((conv(m),m.group(1),m.group(0),pct))
        s=rx.sub(lambda m:' '*len(m.group(0)),s)
    eat(SCI, lambda m: float(m.group(1))*10**int(m.group(2)))
    # A percent sign is a unit: the literal is 100x the underlying fraction, and
    # the registry stores fractions. Accept either reading.
    eat(PCT, lambda m: float(m.group(1)), pct=True)
    eat(DEC, lambda m: float(m.group(1)))
    return out

tot=0; allbad=[]
for f in sys.argv[1:]:
    hits=scan(f); bad=[]
    for v,lit,raw,pct in hits:
        if v==0: continue          # a printed zero has nothing to verify against
        ok = match(v, sig_of(lit)) or (pct and match(v/100.0, sig_of(lit)))
        if not ok: bad.append(raw.strip())
    tot+=len(hits); u=sorted(set(bad))
    allbad+=[(os.path.basename(f),b) for b in u]
    print(f"{os.path.basename(f):28s} {len(hits):4d} literals  {len(bad):3d} unmatched ({len(u)} distinct)")
print(f"\nTOTAL {tot} literals, {len(allbad)} distinct unmatched\n")
for f,b in allbad: print(f"  {f:28s} {b}")
