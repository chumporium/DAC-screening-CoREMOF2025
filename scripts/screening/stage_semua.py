"""Prepare the 7,820 representative CIFs (one per unique material) for the CO2 K_H runs.

- CSD-modified part: extracted from the local zip archive.
- SI part: only the list of names is written; the CIFs are copied from a separate machine.
"""
import csv, io, re, zipfile, os, collections, tarfile

SP = os.path.dirname(os.path.abspath(__file__))
ZIP = os.environ.get('COREMOF_ZIP', 'CSD-modified.zip')      # CoRE MOF 2025 CSD-modified archive


def f(x):
    try:
        return float(x)
    except Exception:
        return None


def basis(n):
    n = re.sub(r'_pacman$', '', n)
    return re.sub(r'_(ASR|FSR|[Ii][Oo][Nn])$', '', n)


# ---- collect all CR rows from both subsets; keep the representative with the largest PLD ----
wakil = {}

for r in csv.DictReader(open(os.path.join(SP, 'si_cr_geom.csv'), encoding='utf-8')):
    p = f(r['PLD'])
    if p is None:
        continue
    b = basis(r['nama'])
    rec = dict(basis=b, sumber='SI', nama=r['nama'], jalur=None, PLD=p,
               air=r['kelas_air'], nama_umum='')
    if b not in wakil or p > wakil[b]['PLD']:
        wakil[b] = rec

z = zipfile.ZipFile(ZIP)
raw = z.read('CSD-modified/CR_data_CSD_modified_20250227.csv').decode('utf-8', 'ignore')
for r in csv.DictReader(io.StringIO(raw)):
    p = f(r['PLD (\u00c5)'])
    if p is None:
        continue
    b = basis(r['refcode'])
    cid = r['coreid']
    ext = 'ASR' if '[ASR]' in cid else ('FSR' if '[FSR]' in cid else 'Ion')
    rec = dict(basis=b, sumber='CSD', nama=r['refcode'],
               jalur='CSD-modified/cifs/CR/%s/%s.cif' % (ext, cid),
               PLD=p, air=r['KH_Classes'], nama_umum=r['name'])
    if b not in wakil or p > wakil[b]['PLD']:
        wakil[b] = rec

print('material unik:', len(wakil))
print('sebaran sumber:', dict(collections.Counter(v['sumber'] for v in wakil.values())))

# ---- write the manifest ----
with open(os.path.join(SP, 'manifes_7820.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=['basis', 'sumber', 'nama', 'jalur', 'PLD', 'air', 'nama_umum'])
    w.writeheader()
    for v in sorted(wakil.values(), key=lambda x: x['basis']):
        w.writerow(v)

# ---- pack the CSD CIFs into a tar.gz archive ----
csd = [v for v in wakil.values() if v['sumber'] == 'CSD']
tgz = os.path.join(SP, 'cif_csd.tar.gz')
n_ok = 0
with tarfile.open(tgz, 'w:gz') as tf:
    for v in csd:
        try:
            data = z.read(v['jalur'])
        except KeyError:
            print('HILANG:', v['jalur']); continue
        info = tarfile.TarInfo(name='cif_semua/%s.cif' % v['nama'].replace('_pacman', ''))
        info.size = len(data)
        tf.addfile(info, io.BytesIO(data))
        n_ok += 1
print('CIF CSD dibungkus:', n_ok, '->', os.path.basename(tgz),
      '%.1f MB' % (os.path.getsize(tgz) / 1e6))

# ---- list of SI names to be copied from the separate machine ----
si = [v['nama'] for v in wakil.values() if v['sumber'] == 'SI']
open(os.path.join(SP, 'daftar_si.txt'), 'w', encoding='utf-8').write('\n'.join(sorted(si)) + '\n')
print('daftar SI ditulis:', len(si), 'nama')
