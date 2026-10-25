"""Collect mixture loadings from RASPA3 output.

The RASPA3 format differs from RASPA2: instead of 'absolute adsorption: ... mol/kg'
it reports 'Abs. loading average  X +/- Y [mol/kg-framework]'.

The header 'Component N (GAS)' appears hundreds of times in one file (once per
statistics block), so taking the first occurrence is wrong. The last occurrence
per gas is used; it holds the production average.

A further pitfall is matching across component boundaries: each loading value
must be assigned to the nearest component header above it, not by order of
appearance.

usage: python3 panen_campuran.py <run_dir> <out.csv> <manifest.csv> [key]
       key = manifest column that identifies the point (e.g. 'indeks' / 'RH')
"""
import os, re, csv, io, sys

PAT_KOMP = re.compile(r'^Component\s+(\d+)\s+\((\w+)\)\s*$', re.M)
PAT_ABS = re.compile(
    r'Abs\. loading average\s+([-\d.eE+]+)\s+\+/-\s+([-\d.eE+]+)\s+\[mol/kg-framework\]')


def baca_satu(path):
    """-> {gas: (q_mol_kg, error)} from the production average."""
    teks = open(path, encoding='utf-8', errors='replace').read()
    judul = [(m.start(), m.group(2)) for m in PAT_KOMP.finditer(teks)]
    if not judul:
        return {}
    nilai = {}
    for m in PAT_ABS.finditer(teks):
        # nearest component header above this value
        gas = None
        for pos, g in judul:
            if pos < m.start():
                gas = g
            else:
                break
        if gas is None:
            continue
        nilai[gas] = (float(m.group(1)), float(m.group(2)))   # last occurrence wins
    return nilai


def panen(dir_lari, keluar, manifes=None, kunci='indeks'):
    man = {}
    if manifes and os.path.exists(manifes):
        for r in csv.DictReader(open(manifes, encoding='utf-8')):
            man[r.get('tag', '')] = r
    baris, gagal = [], []
    for tag in sorted(os.listdir(dir_lari)):
        d = os.path.join(dir_lari, tag)
        od = os.path.join(d, 'output')
        if not os.path.isdir(od):
            continue
        txt = [x for x in os.listdir(od) if x.endswith('.txt')]
        if not txt:
            gagal.append((tag, 'tanpa keluaran')); continue
        nilai = baca_satu(os.path.join(od, txt[0]))
        if not nilai:
            gagal.append((tag, 'tanpa muatan')); continue
        info = man.get(tag, {})
        nama = info.get('nama') or tag.split('__')[0]
        for gas, (q, e) in sorted(nilai.items()):
            baris.append(dict(
                tag=tag, nama=nama, gas=gas, q_mol_kg='%.6e' % q, galat='%.6e' % e,
                galat_pct=('%.1f' % (100 * e / q)) if q > 0 else '',
                **{k: info.get(k, '') for k in
                   ('indeks', 'lapisan', 'RH_setara', 'RH',
                    'p_N2', 'p_O2', 'p_AR', 'p_CO2', 'p_H2O')}))
    if not baris:
        raise SystemExit('tidak ada yang dipanen dari %s' % dir_lari)
    kol = list(baris[0].keys())
    with io.open(keluar, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=kol)
        w.writeheader()
        for b in baris:
            w.writerow(b)
    print('%d baris dari %d lari -> %s'
          % (len(baris), len({b['tag'] for b in baris}), keluar))
    if gagal:
        print('GAGAL %d:' % len(gagal))
        for g in gagal[:10]:
            print('   ', g)
    return baris


if __name__ == '__main__':
    dl, kl = sys.argv[1], sys.argv[2]
    mf = sys.argv[3] if len(sys.argv) > 3 else None
    panen(os.path.expanduser(dl), os.path.expanduser(kl),
          os.path.expanduser(mf) if mf else None)
