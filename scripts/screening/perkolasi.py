"""CO2 accessibility test: height of the energy barrier for traversing the crystal.

E* = lowest energy threshold at which the region {U <= E*} percolates through the crystal.
E* < 0   -> CO2 can enter without crossing any barrier (accessible).
E* large -> pores are sealed in the rigid framework.

Requires only numpy (no scipy). Percolation is detected in a single pass with an
offset-tracking union-find: voxels are added in order of increasing energy, and as
soon as a cluster touches its own periodic image, the energy of that voxel is E*.

usage: python3 perkolasi.py <file.cif> [target_voxels]
"""
import sys, os, re, math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_widom import UFF

KB = 0.0083145            # kJ/mol/K
RC = 6.0                  # steric cutoff; LJ at 6 A is already ~0.01 kJ/mol
CO2 = dict(sigC=2.80, epsC=27.0, sigO=3.05, epsO=79.0, l=1.16)   # TraPPE


def baca(cif):
    txt = open(cif, encoding='utf-8', errors='ignore').read()

    def num(t):
        return float(re.search(r'^\s*' + t + r'\s+([\d.eE+-]+)', txt, re.M).group(1))

    a, b, c = num('_cell_length_a'), num('_cell_length_b'), num('_cell_length_c')
    al, be, ga = [math.radians(num('_cell_angle_' + s)) for s in ('alpha', 'beta', 'gamma')]
    els, fr = [], []
    for l in txt.splitlines():
        s = l.split()
        if len(s) >= 7 and re.match(r'^[A-Z][a-z]?$', s[0]) and re.match(r'^[A-Z]', s[1]):
            try:
                els.append(s[0]); fr.append([float(s[3]), float(s[4]), float(s[5])])
            except ValueError:
                pass
    v1 = [a, 0, 0]
    v2 = [b * math.cos(ga), b * math.sin(ga), 0]
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 1e-12))
    return els, np.array(fr) % 1.0, np.array([v1, v2, [cx, cy, cz]])


def lebar(M):
    V = abs(np.linalg.det(M))
    return np.array([V / np.linalg.norm(np.cross(M[1], M[2])),
                     V / np.linalg.norm(np.cross(M[0], M[2])),
                     V / np.linalg.norm(np.cross(M[0], M[1]))])


def arah(n=13):
    v = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1],
                  [1, 0, -1], [0, 1, 1], [0, 1, -1], [1, 1, 1], [1, 1, -1], [1, -1, 1],
                  [-1, 1, 1]], float)[:n]
    return v / np.linalg.norm(v, axis=1)[:, None]


def peta_energi(cif, target_voxel=60000):
    els, fr, M = baca(cif)
    if not els:
        raise ValueError('CIF tanpa atom')
    hilang = sorted(set(els) - set(UFF))
    if hilang:
        raise ValueError('tanpa UFF: ' + ','.join(hilang))
    L = lebar(M)
    V = abs(np.linalg.det(M))
    kisi = max(0.28, (V / target_voxel) ** (1 / 3.0))
    n = np.maximum(np.ceil(L / kisi).astype(int), 6)

    # only as many periodic atom images as needed
    rep = [int(math.ceil(RC / w)) for w in L]
    off = np.array([[i, j, k] for i in range(-rep[0], rep[0] + 1)
                    for j in range(-rep[1], rep[1] + 1)
                    for k in range(-rep[2], rep[2] + 1)], float)
    A = ((fr[None, :, :] + off[:, None, :]).reshape(-1, 3)) @ M
    epsF = np.tile(np.array([UFF[e][0] for e in els]) * KB, len(off))
    sigF = np.tile(np.array([UFF[e][1] for e in els]), len(off))

    gf = np.stack(np.meshgrid(*[(np.arange(k) + 0.5) / k for k in n], indexing='ij'), -1)
    G = gf.reshape(-1, 3) @ M

    # Cartesian binning with bin size RC
    lo = G.min(axis=0) - RC
    nb = np.maximum(np.ceil((G.max(axis=0) + RC - lo) / RC).astype(int), 1)
    def bid(P):
        b = np.floor((P - lo) / RC).astype(np.int64)
        np.clip(b, 0, nb - 1, out=b)
        return (b[:, 0] * nb[1] + b[:, 1]) * nb[2] + b[:, 2]
    ab = bid(A)
    urut = np.argsort(ab, kind='stable')
    ab_s, A_s, eps_s, sig_s = ab[urut], A[urut], epsF[urut], sigF[urut]
    batas = np.searchsorted(ab_s, np.arange(nb.prod() + 1))

    gb = bid(G)
    situs = [(CO2['sigC'], CO2['epsC'] * KB, 0.0),
             (CO2['sigO'], CO2['epsO'] * KB, +CO2['l']),
             (CO2['sigO'], CO2['epsO'] * KB, -CO2['l'])]
    U = np.full(len(G), 1e6, dtype=np.float32)
    dirs = arah()

    for b in np.unique(gb):
        vid = np.flatnonzero(gb == b)
        bz = b % nb[2]; by = (b // nb[2]) % nb[1]; bx = b // (nb[2] * nb[1])
        cand = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    x, y, zz = bx + dx, by + dy, bz + dz
                    if not (0 <= x < nb[0] and 0 <= y < nb[1] and 0 <= zz < nb[2]):
                        continue
                    q = (x * nb[1] + y) * nb[2] + zz
                    cand.append(np.arange(batas[q], batas[q + 1]))
        if not cand:
            continue
        ci = np.concatenate(cand)
        if ci.size == 0:
            continue
        Ac, ec, sc = A_s[ci], eps_s[ci], sig_s[ci]
        Gv = G[vid]
        terbaik = np.full(len(vid), 1e6, dtype=np.float64)
        for u in dirs:
            tot = np.zeros(len(vid))
            for (sg, ep, d) in situs:
                P = Gv + u * d
                sij = 0.5 * (sc + sg)
                eij = np.sqrt(ec * ep)
                for i in range(0, len(P), 512):
                    dd = np.linalg.norm(P[i:i + 512, None, :] - Ac[None, :, :], axis=2)
                    np.clip(dd, 0.4, None, out=dd)
                    sr6 = (sij[None, :] / dd) ** 6
                    lj = 4 * eij[None, :] * (sr6 * sr6 - sr6)
                    tot[i:i + 512] += np.where(dd < RC, lj, 0.0).sum(axis=1)
            np.minimum(terbaik, tot, out=terbaik)
        U[vid] = terbaik.astype(np.float32)
    return U.reshape(n), n


class UF:
    """Union-find with lattice offsets, used to detect clusters that wrap around the cell."""
    def __init__(self, n):
        self.p = np.arange(n)
        self.off = np.zeros((n, 3), np.int32)

    def cari(self, x):
        o = np.zeros(3, np.int32)
        while self.p[x] != x:
            o = o + self.off[x]
            x = self.p[x]
        return x, o

    def gabung(self, a, b, delta):
        ra, oa = self.cari(a)
        rb, ob = self.cari(b)
        if ra != rb:
            self.p[rb] = ra
            self.off[rb] = oa + delta - ob
            return False
        return not np.array_equal(oa + delta, ob)   # True = wraps around


def ambang(U):
    n1, n2, n3 = U.shape
    N = n1 * n2 * n3
    flat = U.ravel()
    urut = np.argsort(flat, kind='stable')
    ada = np.zeros(N, bool)
    uf = UF(N)
    n23 = n2 * n3
    for v in urut:
        i = v // n23; sisa = v % n23; j = sisa // n3; k = sisa % n3
        ada[v] = True
        for ax, (di, dj, dk) in enumerate([(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                           (0, -1, 0), (0, 0, 1), (0, 0, -1)]):
            ii, jj, kk = i + di, j + dj, k + dk
            d = np.zeros(3, np.int32)
            if ii < 0: ii += n1; d[0] = -1
            elif ii >= n1: ii -= n1; d[0] = +1
            if jj < 0: jj += n2; d[1] = -1
            elif jj >= n2: jj -= n2; d[1] = +1
            if kk < 0: kk += n3; d[2] = -1
            elif kk >= n3: kk -= n3; d[2] = +1
            w = (ii * n2 + jj) * n3 + kk
            if ada[w] and uf.gabung(v, w, d):
                return float(flat[v])
    return None


if __name__ == '__main__':
    cif = sys.argv[1]
    tv = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
    nm = os.path.basename(cif)[:-4]
    try:
        U, n = peta_energi(cif, tv)
        E = ambang(U)
        print('%s\t%.3f\t%s\t%s' % (nm, float(U.min()),
                                    ('%.3f' % E) if E is not None else 'TIDAK_TEMBUS',
                                    'x'.join(map(str, n))))
    except Exception as e:
        print('%s\tGAGAL\t%s\t-' % (nm, str(e)[:60]))
