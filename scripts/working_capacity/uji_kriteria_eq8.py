"""Test of the reachability criterion, Eq. (8), on all 360 working-capacity conditions.

p_des : CO2 partial pressure of the final (self-consistent) desorption gas, y_CO2 x P_des
bound : p_ads x exp[(|dH_CO2|/R)(1/T_ads - 1/T_des)], |dH_CO2| from the adsorption state at 298 K
python uji_kriteria_eq8.py  ->  data/uji_kriteria_eq8.csv
"""
import glob, os
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(ROOT, 'data')


def klas():
    """kerja_final + classification: conv (dq > 2 sd), none (no working capacity), nc (not converged)"""
    a = pd.read_csv(os.path.join(D, 'kerja_final.csv'))
    L = pd.read_csv(os.path.join(D, 'kerja_label_tanpa_kapasitas.csv'))
    L['kls'] = np.where(L.label.str.startswith('tanpa'), 'none', 'nc')
    a = a.merge(L[['nama', 'Ta', 'Td', 'P', 'kls']], on=['nama', 'Ta', 'Td', 'P'], how='left')
    a['kls'] = np.where(a.stabil & a.kls.isna(), np.where(a.dq > 2 * a.sd, 'conv', 'none'), a.kls)
    return a


a = klas()
k = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(D, 'compositions', 'komp_*.csv'))]).drop_duplicates('run').set_index('run')
P = {'p01': 1000.0, 'p05': 5000.0, 'p20': 20265.0}
a['p_des'] = [k.loc[r, 'yCO2'] * P[p] for r, p in zip(a.run, a.P)]
dH = pd.read_csv(os.path.join(D, 'entalpi_ads298.csv')).query('komp == "CO2"').set_index('nama').dH_kJmol.abs()
ads = k[k.index.str.contains('__ads')]


def p_ads(m, Ta):
    rh = 'rh050' if m.endswith('_ion') else 'rh000'
    if Ta == 298:
        return 4.15e-4 * (101325.0 - (1585.0 if rh == 'rh050' else 0.0))
    return ads.loc[f'{m}__ads{Ta}_{rh}', 'yCO2'] * 101325.0


a['p_ads'] = [p_ads(m, Ta) for m, Ta in zip(a.nama, a.Ta)]
a['bound'] = a.p_ads * np.exp(dH.reindex(a.nama).values * 1000 / 8.314 * (1 / a.Ta - 1 / a.Td))
a['below'] = a.p_des < a.bound
a.to_csv(os.path.join(D, 'uji_kriteria_eq8.csv'), index=False)
print(pd.crosstab(a.kls, a.below))
