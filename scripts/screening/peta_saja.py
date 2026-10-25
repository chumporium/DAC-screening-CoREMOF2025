"""Compute only the CO2 energy map and save it as .npy (requires numpy, not scipy).

The percolation test is run separately with scipy (uji_perkolasi.py).
usage: python3 peta_saja.py <file.cif> <dir_out> [target_voxels]
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perkolasi import peta_energi

cif, keluar = sys.argv[1], sys.argv[2]
tv = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
nm = os.path.basename(cif)[:-4]
os.makedirs(keluar, exist_ok=True)
tujuan = os.path.join(keluar, nm + '.npy')
if os.path.exists(tujuan):
    print('%s\tSUDAH_ADA' % nm)
    sys.exit(0)
try:
    U, n = peta_energi(cif, tv)
    np.save(tujuan, U.astype(np.float32))
    print('%s\tOK\t%.3f\t%s' % (nm, float(U.min()), 'x'.join(map(str, n))))
except Exception as e:
    print('%s\tGAGAL\t%s' % (nm, str(e)[:70]))
