"""Henry-based estimate vs five-component GCMC in dry air (Fig. 4, Section 3.4).

Groups: the 40 materials of the GCMC set, the 12 highest-uptake materials that failed the
N2 threshold (Test B), and 18 further rejected materials (8 of Test A, 10 of Test B).
The Henry estimate is K_H,CO2 x p_CO2 at the CO2 partial pressure of the GCMC runs
(415 ppm of 101,325 Pa = 42.05 Pa), so that both refer to the same state.

python korelasi_kh_gcmc.py  ->  data/korelasi_KH_GCMC.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

D = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
b = lambda f: pd.read_csv(os.path.join(D, f))
P_CO2 = 4.15e-4 * 101325.0

kh = b('hasil_screening_co2.csv').set_index('nama')
gc = b('gcmcriil_final.csv').query('rh == "rh000"').set_index('nama').CO2
u = b('ujidac_hasil.csv')
tst = u[(u.uji.isin(['A', 'B'])) & (u.gas == 'CO2') & (u.rh == 'rh000')].set_index('nama').q_mol_kg
b22 = b('uji_B_22.csv').set_index('nama')
top12 = set(b22.index[b22.kelompok == 'serapan_tertinggi'])

rows = []
for n, q in gc.items():
    rows.append((n, 'GCMC set', q))
for n, q in tst.items():
    rows.append((n, 'near N2 threshold (12)' if n in top12 else 'other rejected', q))
k = pd.DataFrame(rows, columns=['nama', 'grup', 'CO2'])
k['KH_CO2'] = kh.KH_CO2.reindex(k.nama).values
k['sigma_KH_pct'] = kh.sigma_pct.reindex(k.nama).values
k['KHpred'] = k.KH_CO2 * P_CO2
k['konvergen'] = k.sigma_KH_pct <= 10
k['rasio'] = k.CO2 / k.KHpred
k.to_csv(os.path.join(D, 'korelasi_KH_GCMC.csv'), index=False)

v = k[k.konvergen & (k.KHpred > 0) & (k.CO2 > 0)]
print(v.grup.value_counts().to_dict(), 'excluded (sigma > 10% or missing):', sorted(set(k.nama) - set(v.nama)))
print(f'n = {len(v)}, Spearman rho = {spearmanr(v.KHpred, v.CO2)[0]:.3f}')
h = v[v.CO2 > 0.1]; print(f'GCMC > 0.1: n = {len(h)}, rho = {spearmanr(h.KHpred, h.CO2)[0]:.3f}')
lo = v[v.KHpred < 0.3]; dev = (lo.rasio - 1).abs()
print(f'Henry estimate < 0.3: {(dev <= 0.2).sum()} of {len(lo)} within 20%, median |dev| {dev.median() * 100:.0f}%')
hi = v[v.KHpred >= 0.3]; print(f'Henry estimate >= 0.3: n = {len(hi)}, ratio {hi.rasio.min():.2f}-{hi.rasio.max():.2f}')
print(v[v.rasio > 1.2][['nama', 'grup', 'KHpred', 'CO2', 'rasio']].to_string())
