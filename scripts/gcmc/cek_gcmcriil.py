import sys
import glob, os, re, statistics as st
T=2.7764
out=[]
for d in sorted(glob.glob(os.path.expanduser(sys.argv[1]+"/*__rh*"))):
    fs=glob.glob(d+'/output/*.txt')
    if not fs: continue
    L=open(fs[0]).read().splitlines()
    if not any('Abs. loading average' in l for l in L): continue
    # production trajectory: instantaneous number of molecules at each print checkpoint
    traj={}; fase=None; comp=None
    for l in L:
        m=re.match(r'(Initialization:|Equilibration:|)s*Current cycle',l)
        if m: fase=(m.group(1) or 'Production').rstrip(':'); continue
        m=re.match(r'Component (\d) \((\w+)\)',l)
        if m: comp=m.group(2); continue
        if fase=='Production' and comp and 'absolute adsorption:' in l and 'molecules' in l:
            traj.setdefault(comp,[]).append(float(l.split()[2])); comp=None
    # final average + blocks (units: molecules/cell)
    names=re.findall(r'^Component (\d) \[(\w+)\]',"\n".join(L),re.M)
    idx=[i for i,l in enumerate(L) if 'Abs. loading average' in l and 'molecules/cell' in l]
    idx=idx[-len(names):]
    for (ci,nm),i in zip(names,idx):
        v=L[i].split(); mean=float(v[3]); pm=float(v[5])
        b=[float(x.split()[-1]) for x in L[i-8:i] if "Block[" in x][-5:]
        sig=pm/T/mean*100 if mean>0 else float('nan')
        med=st.median(b); rasio=max(b)/med if med>0 else float('nan')
        dif=[b[k+1]-b[k] for k in range(4)]
        mono='naik' if all(x>0 for x in dif) else ('turun' if all(x<0 for x in dif) else '-')
        t=traj.get(nm,[]); h=len(t)//2
        paruh=((st.mean(t[h:])-st.mean(t[:h]))/mean*100) if (mean>0 and h>0) else float('nan')
        out.append(f"{os.path.basename(d)},{nm},{mean:.4g},{sig:.2f},{rasio:.2f},{mono},{(b[4]-b[0])/mean*100 if mean>0 else float('nan'):.1f},{paruh:.1f}")
open(os.path.expanduser(sys.argv[2]),'w').write("run,komp,mol_sel,sigma_pct,maks_median,monoton,blok4-blok0_pct,paruh2-paruh1_pct\n"+"\n".join(out)+"\n")
print(len(out))
