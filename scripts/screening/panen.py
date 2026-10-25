"""Collect the CO2 K_H results from ~/khco2/prod and merge them with the manifest."""
import os, re, csv, glob, math

PROD = os.path.expanduser('~/khco2/prod')
MAN = os.path.expanduser('~/manifes_7820.csv')
OUT = os.path.expanduser('~/hasil_kh_co2.csv')

man = {}
for r in csv.DictReader(open(MAN, encoding='utf-8')):
    kunci = r['nama'].replace('_pacman', '')
    man[kunci] = r

pat_kh = re.compile(r'Average Henry coefficient:\s+([\d.eE+-]+|inf|nan)\s+\+/-\s+([\d.eE+-]+|inf|nan)\s+\[mol/kg/Pa\]')
pat_rw = re.compile(r'Average Rosenbluth weight:\s+([\d.eE+-]+|inf|nan)')

baris, gagal = [], []
for d in sorted(os.listdir(PROD)):
    p = os.path.join(PROD, d)
    if not os.path.isdir(p):
        continue
    txts = glob.glob(os.path.join(p, 'output', '*.txt'))
    if not txts:
        gagal.append((d, 'tanpa keluaran')); continue
    isi = open(txts[0], encoding='utf-8', errors='ignore').read()
    m = pat_kh.search(isi)
    if not m:
        gagal.append((d, 'tanpa K_H')); continue
    try:
        kh, sd = float(m.group(1)), float(m.group(2))
    except ValueError:
        gagal.append((d, 'K_H = %s' % m.group(1))); continue
    if not math.isfinite(kh):
        gagal.append((d, 'K_H tak hingga')); continue
    rw = pat_rw.search(isi)
    t = None
    wf = os.path.join(p, 'waktu.txt')
    if os.path.exists(wf):
        for l in open(wf):
            if 'WAKTU_DETIK' in l:
                try:
                    t = float(l.split()[1])
                except Exception:
                    pass
    info = man.get(d, {})
    baris.append(dict(
        nama=d, KH_CO2=kh, sd=sd, galat_pct=100 * sd / kh if kh else '',
        rosenbluth=rw.group(1) if rw else '', detik=t,
        sumber=info.get('sumber', ''), PLD=info.get('PLD', ''),
        air=info.get('air', ''), nama_umum=info.get('nama_umum', '')))

baris.sort(key=lambda x: -x['KH_CO2'])
with open(OUT, 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=['peringkat', 'nama', 'nama_umum', 'KH_CO2', 'sd',
                                       'galat_pct', 'rosenbluth', 'air', 'PLD', 'sumber', 'detik'])
    w.writeheader()
    for i, b in enumerate(baris, 1):
        b['peringkat'] = i
        w.writerow(b)

print('berhasil dipanen :', len(baris))
print('gagal            :', len(gagal))
for g in gagal[:10]:
    print('   ', g)
if baris:
    kh = [b['KH_CO2'] for b in baris]
    q = lambda p: sorted(kh)[int(p * (len(kh) - 1))]
    print()
    print('K_H CO2 (mol/kg/Pa): min %.2e  p25 %.2e  median %.2e  p75 %.2e  p90 %.2e  maks %.2e'
          % (min(kh), q(.25), q(.5), q(.75), q(.9), max(kh)))
    g = [b['galat_pct'] for b in baris if isinstance(b['galat_pct'], float)]
    for amb in (10, 20, 30, 50):
        print('  galat > %2d%% : %5d struktur (%.1f%%)' % (amb, sum(1 for x in g if x > amb),
                                                           100 * sum(1 for x in g if x > amb) / len(g)))
    t = [b['detik'] for b in baris if b['detik']]
    print('  waktu total  : %.1f core-jam, rata-rata %.1f detik/struktur' % (sum(t) / 3600, sum(t) / len(t)))
print('->', OUT)
