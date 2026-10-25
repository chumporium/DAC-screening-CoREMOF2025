"""Percolation test on CO2 energy maps (.npy), using the validated scipy route.

E* = lowest energy threshold at which the region {U <= E*} percolates through the crystal.
  E* < 0            -> accessible; CO2 does not need to cross any barrier
  0 <= E* <~ 25     -> accessible over a thermal barrier (RT at 298 K = 2.48 kJ/mol)
  no percolation    -> pores sealed in the rigid framework
"""
import os, sys, glob, csv
import numpy as np
from scipy import ndimage

STRUK = np.ones((3, 3, 3))     # 26-neighbour connectivity


def menembus(U, t):
    m = U <= t
    if not m.any():
        return False
    T = np.tile(m, (3, 3, 3))
    lab, _ = ndimage.label(T, structure=STRUK)
    n1, n2, n3 = m.shape
    pusat = lab[n1:2 * n1, n2:2 * n2, n3:2 * n3]
    for gs in ((n1, 0, 0), (0, n2, 0), (0, 0, n3)):
        g = lab[n1 + gs[0]:2 * n1 + gs[0], n2 + gs[1]:2 * n2 + gs[1], n3 + gs[2]:2 * n3 + gs[2]]
        if ((pusat > 0) & (pusat == g)).any():
            return True
    return False


def ambang(U, batas=200.0):
    if not menembus(U, batas):
        return None
    lo, hi = float(U.min()), batas
    if menembus(U, lo):
        return lo
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if menembus(U, mid):
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == '__main__':
    src = sys.argv[1]
    out = sys.argv[2]
    hasil = []
    berkas = sorted(glob.glob(os.path.join(src, '*.npy')))
    print('peta ditemukan:', len(berkas))
    for i, f in enumerate(berkas, 1):
        nm = os.path.basename(f)[:-4]
        U = np.load(f)
        E = ambang(U)
        hasil.append(dict(nama=nm, U_min=round(float(U.min()), 3),
                          E_halangan=('%.3f' % E) if E is not None else '',
                          terjangkau=('tersegel' if E is None else
                                      ('bebas' if E < 0 else ('termal' if E < 25 else 'berat'))),
                          voxel='x'.join(map(str, U.shape))))
        if i % 50 == 0:
            print('  %d/%d' % (i, len(berkas)), flush=True)
    with open(out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['nama', 'U_min', 'E_halangan', 'terjangkau', 'voxel'])
        w.writeheader()
        for h in hasil:
            w.writerow(h)
    import collections
    print('ringkasan:', dict(collections.Counter(h['terjangkau'] for h in hasil)))
    print('->', out)
