"""Stage-1 screening: pore size accessible to CO2 (kinetic diameter 3.30 A).

Merges the two CR subsets of CoRE MOF 2025 (SI + CSD-modified), removes
ASR/FSR/Ion duplicates, then counts the structures passing several PLD thresholds.
"""
import csv, io, re, zipfile, collections, math, json, sys, os

SP = os.environ.get('WORK_DIR', os.path.dirname(os.path.abspath(__file__)))      # working directory
ZIP = os.environ.get('COREMOF_ZIP', 'CSD-modified.zip')      # CoRE MOF 2025 CSD-modified archive

def f(x):
    try:
        return float(x)
    except Exception:
        return None

baris = []

# --- SI subset ---
for r in csv.DictReader(open(os.path.join(SP, 'si_cr_geom.csv'), encoding='utf-8')):
    baris.append({'sumber': 'SI', 'nama': r['nama'], 'PLD': f(r['PLD']), 'LCD': f(r['LCD']),
                  'rho': f(r['density']), 'VF': f(r['VF']), 'PV': f(r['PV_cm3g']),
                  'ASA': f(r['ASA_m2g']), 'air': r['kelas_air'], 'KH_air': f(r['KH_air']),
                  'natom': f(r['n_atoms']), 'logam': r['logam'], 'nama_umum': ''})

# --- CSD-modified subset ---
z = zipfile.ZipFile(ZIP)
raw = z.read('CSD-modified/CR_data_CSD_modified_20250227.csv').decode('utf-8', 'ignore')
for r in csv.DictReader(io.StringIO(raw)):
    baris.append({'sumber': 'CSD', 'nama': r['refcode'], 'PLD': f(r['PLD (\u00c5)']),
                  'LCD': f(r['LCD (\u00c5)']), 'rho': f(r['Density (g/cm3)']), 'VF': f(r['VF']),
                  'PV': f(r['PV (cm3/g)']), 'ASA': f(r['ASA (m2/g)']), 'air': r['KH_Classes'],
                  'KH_air': None, 'natom': f(r['natoms']), 'logam': r['Metal Types'],
                  'nama_umum': r['name']})

print('baris CR total          :', len(baris))
print('  SI                    :', sum(1 for b in baris if b['sumber'] == 'SI'))
print('  CSD-modified          :', sum(1 for b in baris if b['sumber'] == 'CSD'))
tanpa = [b for b in baris if b['PLD'] is None]
print('tanpa nilai PLD         :', len(tanpa))

# --- unique material basis: strip the ASR/FSR/Ion and _pacman suffixes ---
def basis(n):
    n = re.sub(r'_pacman$', '', n)
    n = re.sub(r'_(ASR|FSR|[Ii][Oo][Nn])$', '', n)
    return n

grup = collections.defaultdict(list)
for b in baris:
    if b['PLD'] is None:
        continue
    b['basis'] = basis(b['nama'])
    grup[b['basis']].append(b)
print('material unik (buang varian ASR/FSR/Ion):', len(grup))

# representative of each material: the variant with the largest PLD (most open after curation)
wakil = []
for k, v in grup.items():
    wakil.append(max(v, key=lambda x: x['PLD']))
print('wakil dipakai           :', len(wakil))
print()

# --- PLD distribution ---
plds = sorted(b['PLD'] for b in wakil)
def kuantil(p):
    return plds[min(len(plds) - 1, int(p * len(plds)))]
print('sebaran PLD (A): min %.2f  p10 %.2f  p25 %.2f  median %.2f  p75 %.2f  p90 %.2f  maks %.2f'
      % (plds[0], kuantil(.10), kuantil(.25), kuantil(.50), kuantil(.75), kuantil(.90), plds[-1]))
print()

# --- thresholds ---
AMBANG = [(2.4, 'ambang porositas CoRE'),
          (3.3, 'diameter kinetik CO2 (3,30 A) - batas mutlak'),
          (3.8, 'CO2 + 0,5 A - difusi wajar'),
          (4.3, 'CO2 + 1,0 A - longgar'),
          (5.0, 'pori sedang')]
print('%-8s %-38s %8s %8s' % ('PLD >=', 'arti', 'lolos', '%'))
for a, ket in AMBANG:
    n = sum(1 for b in wakil if b['PLD'] >= a)
    print('%-8.1f %-38s %8d %7.1f%%' % (a, ket, n, 100 * n / len(wakil)))
print()

# --- cross-tabulation with hydrophobicity ---
HID = ('weak', 'superweak')
print('%-8s %10s %12s %14s' % ('PLD >=', 'lolos', 'hidrofobik', '% hidrofobik'))
for a, _ in AMBANG:
    sub = [b for b in wakil if b['PLD'] >= a]
    h = [b for b in sub if b['air'] in HID]
    print('%-8.1f %10d %12d %13.1f%%' % (a, len(sub), len(h), 100 * len(h) / max(1, len(sub))))
print()

# --- zero ASA: pores inaccessible to the N2 probe even though the PLD passes ---
lolos = [b for b in wakil if b['PLD'] >= 3.3]
asa0 = [b for b in lolos if b['ASA'] is not None and b['ASA'] == 0]
print('lolos PLD>=3,3 tapi ASA = 0 m2/g :', len(asa0), '(%.1f%%)' % (100 * len(asa0) / len(lolos)))

# --- working set: PLD >= 3.8 A and hydrophobic ---
kerja = [b for b in wakil if b['PLD'] >= 3.8 and b['air'] in HID]
print()
print('HIMPUNAN KERJA (PLD >= 3,8 A dan hidrofobik):', len(kerja))
print('  sumber:', dict(collections.Counter(b['sumber'] for b in kerja)))
print('  logam terbanyak:', collections.Counter(b['logam'] for b in kerja).most_common(8))
bern = [b for b in kerja if b['nama_umum'] not in ('', '-', None)]
print('  punya nama umum:', len(bern), [b['nama_umum'] for b in bern[:12]])

with open(os.path.join(SP, 'himpunan_kerja.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=list(kerja[0].keys()))
    w.writeheader()
    for b in sorted(kerja, key=lambda x: -x['PLD']):
        w.writerow(b)
print('  -> ditulis ke himpunan_kerja.csv')
