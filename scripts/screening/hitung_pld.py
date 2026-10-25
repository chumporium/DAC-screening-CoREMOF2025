"""Independent PLD calculation: hard spheres + percolation test, with interchangeable radius tables.

PLD = diameter of the largest sphere that can traverse the crystal (pass from one face to
the opposite face through periodic images). Computed on a voxel grid with connected-component
labelling on a 3x3x3 supercell, then checking whether a component touches itself
after a shift by one lattice vector.
"""
import sys, re, math
import numpy as np
from scipy import ndimage

# --- radius tables (Angstrom) ---
RADII = {
    'bondi': {'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'Si': 2.10, 'S': 1.80,
              'Cl': 1.75, 'Zn': 1.39, 'Cu': 1.40, 'Ni': 1.63, 'Nb': 2.00, 'Fe': 2.00,
              'Co': 2.00, 'Cd': 1.58, 'Ag': 1.72, 'V': 2.00, 'Sc': 2.00},
    # UFF: x_i (vdW distance, r_min) divided by 2
    'uff':   {'H': 1.443, 'C': 1.926, 'N': 1.830, 'O': 1.750, 'F': 1.682, 'Si': 2.147,
              'S': 2.017, 'Cl': 1.973, 'Zn': 1.382, 'Cu': 1.748, 'Ni': 1.417, 'Nb': 1.582,
              'Fe': 1.456, 'Co': 1.436, 'Cd': 1.424, 'Ag': 1.574, 'V': 1.572, 'Sc': 1.648},
    # Cordero (covalent radii) - much smaller
    'kovalen': {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Si': 1.11, 'S': 1.05,
                'Cl': 1.02, 'Zn': 1.22, 'Cu': 1.32, 'Ni': 1.24, 'Nb': 1.64, 'Fe': 1.32,
                'Co': 1.26, 'Cd': 1.44, 'Ag': 1.45, 'V': 1.53, 'Sc': 1.70},
    # Zeo++ default table (radii.rad) - vdW, but H smaller than Bondi
    'zeopp': {'H': 1.09, 'C': 1.67, 'N': 1.60, 'O': 1.52, 'F': 1.47, 'Si': 2.10, 'S': 1.80,
              'Cl': 1.75, 'Zn': 1.39, 'Cu': 1.40, 'Ni': 1.63, 'Nb': 2.00, 'Fe': 2.00,
              'Co': 2.00, 'Cd': 1.58, 'Ag': 1.72, 'V': 2.00, 'Sc': 2.00},
}


def baca(cif):
    txt = open(cif, encoding='utf-8', errors='ignore').read()

    def num(t):
        return float(re.search(r'^\s*' + t + r'\s+([\d.eE+-]+)', txt, re.M).group(1))

    a, b, c = num('_cell_length_a'), num('_cell_length_b'), num('_cell_length_c')
    al, be, ga = [math.radians(num('_cell_angle_' + s)) for s in ('alpha', 'beta', 'gamma')]
    els, fr = [], []
    for l in txt.splitlines():
        s = l.split()
        if len(s) == 8 and re.match(r'^[A-Z][a-z]?$', s[0]) and re.match(r'^[A-Z][a-z]?\d', s[1]):
            els.append(s[0]); fr.append([float(s[3]), float(s[4]), float(s[5])])
    # cell matrix (rows = lattice vectors)
    v1 = [a, 0, 0]
    v2 = [b * math.cos(ga), b * math.sin(ga), 0]
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 1e-12))
    M = np.array([v1, v2, [cx, cy, cz]])
    return els, np.array(fr) % 1.0, M


def peta_ruang(cif, tabel, kisi=0.20):
    """Return the array r_free[i,j,k] = radius of the largest sphere that fits at that voxel."""
    els, fr, M = baca(cif)
    R = np.array([RADII[tabel].get(e, 1.7) for e in els])
    L = np.linalg.norm(M, axis=1)
    n = np.maximum(np.ceil(L / kisi).astype(int), 8)
    gf = np.stack(np.meshgrid(*[(np.arange(k) + 0.5) / k for k in n], indexing='ij'), -1)
    G = gf.reshape(-1, 3) @ M  # Cartesian voxel coordinates
    # periodic images of the atoms: 3x3x3
    off = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    Af = (fr[None, :, :] + off[:, None, :]).reshape(-1, 3)
    A = Af @ M
    Rall = np.tile(R, len(off))
    best = np.full(len(G), 1e9)
    blok = 4000
    for i in range(0, len(G), blok):
        d = np.linalg.norm(G[i:i + blok, None, :] - A[None, :, :], axis=2) - Rall[None, :]
        best[i:i + blok] = d.min(axis=1)
    return best.reshape(n), n, M, len(els)


def menembus(bebas3d, rp):
    """Does the region with r_free >= rp percolate through the crystal?"""
    m = bebas3d >= rp
    if not m.any():
        return False
    T = np.tile(m, (3, 3, 3))
    lab, _ = ndimage.label(T, structure=np.ones((3, 3, 3)))
    nx, ny, nz = m.shape
    # voxels in the central cell
    pusat = lab[nx:2 * nx, ny:2 * ny, nz:2 * nz]
    for gs in [(nx, 0, 0), (0, ny, 0), (0, 0, nz)]:
        geser = lab[nx + gs[0]:2 * nx + gs[0], ny + gs[1]:2 * ny + gs[1], nz + gs[2]:2 * nz + gs[2]]
        sama = (pusat > 0) & (pusat == geser)
        if sama.any():
            return True
    return False


def pld(cif, tabel, kisi=0.20):
    bebas, n, M, na = peta_ruang(cif, tabel, kisi)
    lo, hi = 0.0, float(bebas.max())
    if not menembus(bebas, 1e-6):
        return 0.0, na, n
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if menembus(bebas, mid):
            lo = mid
        else:
            hi = mid
    return 2 * lo, na, n


if __name__ == '__main__':
    kisi = float(sys.argv[1]) if len(sys.argv) > 1 else 0.20
    import glob, os
    SP = os.path.dirname(os.path.abspath(__file__))
    berkas = sorted(glob.glob(os.path.join(SP, 'uji_co2_asli', '*.cif')))
    print('%-12s %6s %-10s %8s %s' % ('struktur', 'atom', 'tabel', 'PLD', 'kisi voxel'))
    for f in berkas:
        nm = os.path.basename(f)[:-4]
        for tabel in ('uff', 'bondi', 'zeopp', 'kovalen'):
            p, na, n = pld(f, tabel, kisi)
            print('%-12s %6d %-10s %8.2f %s' % (nm, na, tabel, p, 'x'.join(map(str, n))))
