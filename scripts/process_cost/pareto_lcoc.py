"""Two-objective screening of operating conditions: minimum levelised cost of capture (LCOC) and maximum
CO2 purity, per site and sorbent, by exhaustive enumeration of a small discrete design space.

Decision variables (desorption pressure fixed at 1 kPa, where the desorption surrogate is validated):
  T_des                         348, 373, 398 K
  drying target (scenario B)    RH at the contactor 0 (dew point of the dryer), 1, 2, 5, 10 %
Sorbents: KENVOR (scenario A), SIFSIX-3-Ni (scenario B), SIFSIX-3-Ni scaled to the experimental uptake.
The plant model dacsys is used read-only through dacsys_gcmc.py (energy supply sized for each case).

python pareto_lcoc.py [n_proc]  -> data/pareto_lcoc.csv (all cases, column 'pareto' marks non-dominated ones)
"""
import os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import dacsys_gcmc as dg
from dacsys import plant, sorbents, climate_sites as cs
from dacsys.levers import Case, Prm

T_DES = (348, 373, 398)
RH_B = (0.0, 0.01, 0.02, 0.05, 0.10)
BAHAN = {'KENVOR_ion': ('A', None), 'OCEGOV_ASR': ('B', None), 'OCEGOV_ASR:kor': ('B', 0.49)}
for Td in T_DES:                                            # in-memory regeneration settings, no file written
    plant.REGEN[f'gcmc{Td}'] = (Td + 0.15, 0.01, 'vacuum')


def daftar():
    out = []
    for key, (sk, f) in BAHAN.items():
        m = key.split(':')[0]
        for rh in ((1.0,) if sk == 'A' else RH_B):
            name = f'PARETO:{key}:rh{rh:.2f}'
            s = dg.DummySorbent(name, ((1.0, -20.0, 40.0),), hoa=abs(float(dg.h298.loc[m, 'CO2'])), T_des=373.15, p_des=0.01,
                                mode='vacuum', cost=dg.PRODUKSI, hydrophobic=(sk == 'A'), rh_max=rh,
                                cp=float(dg.cp.loc[m].Cp_300K), rho=dg.rho_of(m), meta={'nama': m, 'f': f})
            sorbents.register(s)
            out.append((key, sk, rh, name))
    return out


def jalankan(site):
    w = cs.load(site)
    res = plant.resources(w, dg.PRM)
    rows = []
    for key, sk, rh, name in daftar():
        for Td in T_DES:
            prm = Prm(purity_model='henry', regen=f'gcmc{Td}', p_des=0.01)
            c = Case(site=site, drying='none' if sk == 'A' else 'cond_silica', cold='ambient', sorbent=name, lever='PARETO', prm=prm)
            r = plant.evaluate(c, w, res)
            rows.append(dict(site=site, sorben=key, skenario=sk, T_des=Td, rh_target=rh,
                             **{k: r.get(k) for k in ('status', 'lcoc', 'purity', 'wcap_mmolg', 'hours', 'heat_GJ', 'heat_sil_GJ', 'elec_GJ',
                                                      'capex_sorbent', 'capex_process', 'capex_energy', 'co2_t_per_kgs')}))
    return rows


def pareto(g):
    """non-dominated rows: lower LCOC and higher purity"""
    ok = g[g.status == 'ok'].sort_values(['lcoc', 'purity'], ascending=[True, False])
    best, idx = -np.inf, []
    for i, r in ok.iterrows():
        if r.purity > best + 1e-9:
            idx.append(i); best = r.purity
    return idx


if __name__ == '__main__':
    from multiprocessing import Pool
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    with Pool(n) as pool:
        rows = sum(pool.map(jalankan, list(cs.SITES)), [])
    out = pd.DataFrame(rows)
    tm = pd.read_csv(os.path.join(dg.D, 'lcoc_gcmc.csv')).groupby('site').T_mean_C.first()
    out['T_mean_C'] = out.site.map(tm)
    out['pareto'] = False
    for _, g in out.groupby(['site', 'sorben']):
        out.loc[pareto(g), 'pareto'] = True
    out.to_csv(os.path.join(dg.D, 'pareto_lcoc.csv'), index=False)
    best = out[out.status == 'ok'].sort_values('lcoc').groupby(['site', 'sorben']).head(1)
    pd.set_option('display.width', 220)
    print(best.pivot_table(index='site', columns='sorben', values='lcoc').round(0).to_string())
    print(best.groupby('sorben')[['T_des', 'rh_target']].agg(lambda s: s.value_counts().to_dict()))
