"""Collect the working-capacity runs: uptake (mol/kg) and adsorption enthalpy per component.
usage: python3 panen_kerja.py <out.csv>"""
import glob, os, re, sys
NAN = float("nan"); T = 2.7764
out = ["run,komp,q_molkg,sigma_pct,dH_kJmol,sig_dH"]
for d in sorted(glob.glob(os.path.expanduser("~/kerja/*__*"))):
    fs = glob.glob(d + "/output/*.txt")
    if not fs: continue
    s = open(fs[0]).read()
    if "Abs. loading average" not in s: continue
    L = s.splitlines()
    nm = [m[1] for m in re.findall(r"^Component (\d) \[(\w+)\]", s, re.M)]
    v = [l.split() for l in L if "Abs. loading average" in l and "mol/kg-framework" in l][-len(nm):]
    h = {}
    i0 = s.rfind("Enthalpy of adsorption\n")
    if i0 >= 0:
        for c, a, b in re.findall(r"Component \d \[(\w+)\].*?\n\s+(-?[\d.e+-]+) \+/-\s+([\d.e+-]+) \[kJ/mol\]", s[i0:], re.S):
            h.setdefault(c, (float(a), float(b)))
    for n, x in zip(nm, v):
        q = float(x[3]); pm = float(x[5]); sg = pm / T / q * 100 if q > 0 else NAN
        dh, sdh = h.get(n, (NAN, NAN))
        out.append("%s,%s,%.5g,%.1f,%.4g,%.3g" % (os.path.basename(d), n, q, sg, dh, sdh / T))
open(sys.argv[1], "w").write("\n".join(out) + "\n")
print(len(out) - 1)
