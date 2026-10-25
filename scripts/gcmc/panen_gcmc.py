# Harvest of five-component GCMC runs in ~/<DIR>/<material>__rhXXX (default uji415): final block averages where the run
# finished; for runs stopped early, the running production average (truncated) or, for runs
# still in equilibration, the mean of the instantaneous loadings over the last 10 print blocks.
import re, glob, os, statistics, sys
DIR = sys.argv[1] if len(sys.argv) > 1 else 'uji415'
KOMP = ['N2', 'O2', 'AR', 'H2O', 'CO2']
out = ['run,status,cycles_prod,' + ','.join(f'{k}_molkg,{k}_ci95' for k in KOMP)]
for d in sorted(glob.glob(os.path.expanduser('~/' + DIR + '/*__rh*/'))):
    run = os.path.basename(d.rstrip('/'))
    f = sorted(glob.glob(d + 'output/*.txt'))[0]
    s = open(f, errors='ignore').read()
    val = {}
    if 'Abs. loading average' in s:
        nm = [m[1] for m in re.findall(r'^Component (\d) \[(\w+)\]', s, re.M)]
        L = [l.split() for l in s.splitlines() if 'Abs. loading average' in l and 'mol/kg-framework' in l][-len(nm):]
        for k, l in zip(nm, L):
            val[k] = (float(l[3]), float(l[5]))
        status, ncyc = 'finished', 50000
    else:
        blocks = re.split(r'\n(?=(?:\w+: )?Current cycle: )', s)[1:]
        last = blocks[-1]
        m = re.match(r'(?:(\w+): )?Current cycle: (\d+)', last)
        if m.group(1) is None:
            status, ncyc = 'stopped_in_production_running_average', int(m.group(2))
            for k in KOMP:
                mm = re.search(r'Component \d \(' + k + r'\)\s+absolute adsorption:.*?\n.*?\n\s+[-\d.e+]+ mol/kg\s+\(([-\d.e+]+) \+/- ([-\d.e+]+)\)', last)
                if mm: val[k] = (float(mm.group(1)), float(mm.group(2)))
        else:
            status, ncyc = f'stopped_in_{m.group(1).lower()}_mean_last10_blocks', 0
            for k in KOMP:
                xs = []
                for b in blocks[-10:]:
                    mm = re.search(r'Component \d \(' + k + r'\)\s+absolute adsorption:.*?\n.*?\n\s+([-\d.e+]+) mol/kg', b)
                    if mm: xs.append(float(mm.group(1)))
                if xs: val[k] = (statistics.mean(xs), float('nan'))
    row = [run, status, str(ncyc)]
    for k in KOMP:
        a, b = val.get(k, (float('nan'), float('nan')))
        row += [f'{a:.6e}', f'{b:.6e}']
    out.append(','.join(row))
open(os.path.expanduser('~/' + DIR + '_hasil.csv'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
