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
        if abs(av-aw) <= 0.5*ulp*1.001: out.append(k)
    return out

SCI=re.compile(r'([+-]?\d+\.?\d*)\s*\\times\s*10\^\{(-?\d+)\}')
PCT=re.compile(r'([+-]?\d+\.\d{2,})\s*\\%')
DEC=re.compile(r'(?<![\d.])(\d+\.\d{3,})(?!\d)')

def scan(path):
    s=re.sub(r'(?m)^%.*','',open(path).read())
    # LaTeX lengths and font sizes are typography, not data
    s=re.sub(r'\{[\d.]+\\(?:text|line|column|page)width\}', ' ', s)
    s=re.sub(r'\\(?:setlength|hspace|vspace|tabcolsep|emergencystretch)[^\n]*', ' ', s); out=[]
    def eat(rx,conv):
        nonlocal s
        for m in list(rx.finditer(s)): out.append((conv(m),m.group(1),m.group(0)))
        s=rx.sub(lambda m:' '*len(m.group(0)),s)
    eat(SCI, lambda m: float(m.group(1))*10**int(m.group(2)))
    eat(PCT, lambda m: float(m.group(1)))
    eat(DEC, lambda m: float(m.group(1)))
    return out

tot=0; allbad=[]
for f in sys.argv[1:]:
    hits=scan(f); bad=[]
    for v,lit,raw in hits:
        if v==0: continue          # a printed zero has nothing to verify against
        if not match(v, sig_of(lit)): bad.append(raw.strip())
    tot+=len(hits); u=sorted(set(bad))
    allbad+=[(os.path.basename(f),b) for b in u]
    print(f"{os.path.basename(f):28s} {len(hits):4d} literals  {len(bad):3d} unmatched ({len(u)} distinct)")
print(f"\nTOTAL {tot} literals, {len(allbad)} distinct unmatched\n")
for f,b in allbad: print(f"  {f:28s} {b}")
