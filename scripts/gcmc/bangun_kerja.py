"""Build GCMC inputs from a composition table: run,T,P,yN2,yO2,yAR,yH2O,yCO2
usage: python3 bangun_kerja.py <composition.csv> <dir_out> <run_list.txt>"""
import sys, os, csv, json, glob
H = os.path.expanduser('~'); sys.path.insert(0, H)
from gen_pindai import tulis
DIRS = ['raspa_skrining/cif', 'cif_pjg468', 'raspa_skrining/cif_kolam', 'cif219'] + \
       [os.path.basename(p) for p in glob.glob(H + '/cif_*') if os.path.isdir(p)]
def cari(nm):
    for d in DIRS:
        p = f'{H}/{d}/{nm}.cif'
        if os.path.exists(p): return p
    raise FileNotFoundError(nm)
src, out, daftar = sys.argv[1:4]
L = []
for r in csv.DictReader(open(src)):
    run = r['run']; nm = run.split('__')[0]; dst = f'{out}/{run}'
    y = {g: float(r['y' + g]) for g in ('N2', 'O2', 'AR', 'H2O', 'CO2')}
    tulis(cari(nm), dst, 50, siklus=50000, suhu=float(r['T']), grid=False)
    s = json.load(open(dst + '/simulation.json'))
    s['NumberOfInitializationCycles'] = 50000; s['NumberOfEquilibrationCycles'] = 50000
    s['Systems'][0]['ExternalTemperature'] = float(r['T'])
    s['Systems'][0]['ExternalPressure'] = float(r['P'])
    kom = []
    for c in s['Components']:
        if y.get(c['Name'], 0) > 0:
            c['MolFraction'] = y[c['Name']]; kom.append(c)
        else:
            f = f"{dst}/{c['Name']}.json"
            if os.path.exists(f): os.remove(f)
    s['Components'] = kom
    json.dump(s, open(dst + '/simulation.json', 'w'), indent=1)
    L.append(run)
open(daftar, 'w').write('\n'.join(L) + '\n')
print(len(L), 'lari')
