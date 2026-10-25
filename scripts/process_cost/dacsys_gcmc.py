"""Levelised cost of capture with GCMC-derived sorbent behaviour and the dacsys plant/energy/cost model.

The dacsys package (path in the environment variable DACSYS_PATH) is used read-only: no file in that folder is written
(no bytecode, all output to this project). Its air train supplies the air-side process (drying, fans,
optional chilling); every sorbent-related quantity is replaced by the five-component GCMC map of this
work (data/peta_iklim.csv, desorption at 1 kPa):
  working capacity(T_ads, RH at the contactor), CO2 purity of the released gas, desorption heat of all
  released components, vacuum work for all released gases, void evacuation work, co-released water.
Sensible heat, cycle time (heat transfer), energy block and capital costs are dacsys'.

python dacsys_gcmc.py uji            -> one site, quick check
python dacsys_gcmc.py semua [n_proc] -> 21 sites x cases -> data/lcoc_gcmc.csv
PERLAKUAN=inert|umpan python dacsys_gcmc.py semua -> data/lcoc_gcmc_inert.csv / lcoc_gcmc_umpan.csv (simplified treatments)
"""
import os, sys
sys.dont_write_bytecode = True
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
ARS = os.environ.get('DACSYS_PATH', 'dacsys_path_not_set')      # folder that contains the dacsys package
sys.path.insert(0, ARS)

from dataclasses import replace
import numpy as np, pandas as pd

from dacsys import plant, sorbents, tvsa, climate_sites as cs
from dacsys.levers import Case, Prm
from dacsys.weather import humidity_ratio

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(ROOT, 'data')
M_CO2, M_H2O = 44.01, 18.015
KJMOL_PER_GJT = 1.0 / (1000.0 / M_CO2 / 1000.0)        # GJ/t -> kJ/mol CO2

# ---------------------------------------------------------------- sorbent data (this work)
cp = pd.read_csv(os.path.join(D, 'cp_semua.csv')).set_index('nama')
geo = pd.read_csv(os.path.join(D, 'saring_tahap1.csv'))
geo['nama'] = geo.nama.str.replace('_pacman', '', regex=False)
geo = geo.set_index('nama')
h298 = pd.read_csv(os.path.join(D, 'entalpi_ads298.csv')).pivot_table(index='nama', columns='komp', values='dH_kJmol')
EUR_USD = 1.0813     # 2023 average reference rate (ECB), USD per EUR
PRODUKSI = 30.0 * EUR_USD   # $/kg, uniform material specific cost of the PrISMa platform (30 EUR/kg) for every MOF
UMUR_SORBEN = 5      # years, sorbent replaced every 5 years (PrISMa material life time)


def faktor_ganti(p):
    """Present value of the initial sorbent plus replacements every UMUR_SORBEN years, relative to the initial sorbent."""
    return sum((1.0 + p.discount) ** (-k) for k in range(0, p.lifetime, UMUR_SORBEN))
BAHAN = {  # name: (scenario, tolerated RH fraction at T_ads, experimental/GCMC uptake ratio or None)
    'KENVOR_ion': ('A', 1.0, None), 'KENVOR02_ion': ('A', 1.0, None),
    'OCEGOV_ASR': ('B', 0.05, 0.49), 'ARAHIM01_ASR': ('B', 0.0, 0.73),
}


def rho_of(m):
    return float(geo.loc[m].rho) if m in geo.index else float(geo[geo.basis == m.split('_')[0][:6]].iloc[0].rho)


class DummySorbent(sorbents.Sorbent):
    """Stand-in sorbent that keeps the feasibility flags of dacsys neutral; all sorbent quantities are overwritten."""
    def working_capacity(self, T_ads, p_co2_bar):
        return np.full_like(np.asarray(T_ads, float), 1.0e-3 * M_CO2)

    def selectivity_n2(self, T):
        return np.full_like(np.asarray(T, float), 1e12)


def daftar_sorben(koreksi=False):
    out = {}
    for m, (sk, rh_tol, f) in BAHAN.items():
        name = f'GCMC:{m}' + (':kor' if (koreksi and f) else '')
        s = DummySorbent(name, ((1.0, -20.0, 40.0),), hoa=abs(float(h298.loc[m, 'CO2'])), T_des=373.15, p_des=0.01,
                         mode='vacuum', cost=PRODUKSI, hydrophobic=(sk == 'A'), rh_max=rh_tol,
                         cp=float(cp.loc[m].Cp_300K), rho=rho_of(m), meta={'nama': m, 'f': f if koreksi else None})
        sorbents.register(s)
        out[name] = s
    return out


# ---------------------------------------------------------------- GCMC map interpolation
PERLAKUAN = os.environ.get('PERLAKUAN', 'lima')     # 'inert' / 'umpan': simplified desorption treatments (peta_energi.py)
PETA = pd.read_csv(os.path.join(D, 'peta_iklim.csv' if PERLAKUAN == 'lima' else f'peta_iklim_{PERLAKUAN}.csv'))


class Peta:
    """Quantities per mol CO2 released as functions of (T_ads, RH) at a given T_des."""
    def __init__(self, m, Td=373, f=None):
        x = PETA[(PETA.nama == m) & (PETA.Td == Td) & (PETA.dq_CO2 > 0)].copy()
        if f:
            x['dq_CO2'] = x.dq_CO2 * f
        x['n_rel'] = (x.dq_CO2 + x.H2O_lepas.fillna(0) + x.N2O2Ar_lepas.fillna(0)) / x.dq_CO2
        x['w_rel'] = x.H2O_lepas.fillna(0) / x.dq_CO2
        x['qdes'] = x.panas_desorpsi * KJMOL_PER_GJT * (1.0 if not f else 1.0)   # kJ/mol CO2, all released gases
        self.RH = sorted(x.RH.unique())
        self.tab = {r: x[x.RH == r].sort_values('Ta') for r in self.RH}

    def _kurva(self, r, col, T, log=False):
        t = self.tab[r]
        Tn, v = t.Ta.values.astype(float), t[col].values.astype(float)
        if len(Tn) == 1:
            return np.full_like(T, v[0])
        y = np.log(np.maximum(v, 1e-9)) if log else v
        out = np.interp(1.0 / T, 1.0 / Tn[::-1], y[::-1])
        hi = T > Tn.max()                                   # warmer than grid: extrapolate in 1/T from last two nodes
        if hi.any():
            s = (y[-1] - y[-2]) / (1.0 / Tn[-1] - 1.0 / Tn[-2])
            out = np.where(hi, y[-1] + s * (1.0 / T - 1.0 / Tn[-1]), out)
        return np.exp(out) if log else out

    def __call__(self, col, T, RH, log=False):
        T = np.clip(np.asarray(T, float), 243.0, None)
        RH = np.clip(np.asarray(RH, float) * 100.0, 0, 100)
        rr = np.array(self.RH, float)
        vals = np.stack([self._kurva(r, col, T, log) for r in self.RH])     # (nRH, H)
        idx = np.clip(np.searchsorted(rr, RH) - 1, 0, len(rr) - 2) if len(rr) > 1 else np.zeros_like(RH, int)
        if len(rr) == 1:
            return vals[0]
        r0, r1 = rr[idx], rr[idx + 1]
        wgt = np.clip((RH - r0) / np.maximum(r1 - r0, 1e-9), 0, 1)
        H = np.arange(len(T))
        a, b = vals[idx, H], vals[idx + 1, H]
        return np.exp((1 - wgt) * np.log(np.maximum(a, 1e-9)) + wgt * np.log(np.maximum(b, 1e-9))) if log else (1 - wgt) * a + wgt * b


_ASLI = plant.air_train
_PETA = {}


def air_train_gcmc(w, case, T_design, dT_rec=3.0):
    L = _ASLI(w, case, T_design, dT_rec)
    p, sb = case.prm, plant.sorbent_for(case)
    m, f = sb.meta['nama'], sb.meta.get('f')
    Td = int(round(sb.T_des))                                  # GCMC map at the desorption temperature of the case
    key = (m, f, Td)
    if key not in _PETA:
        _PETA[key] = Peta(m, Td, f)
    pt = _PETA[key]
    P = w['P']
    # after-cooler behind an ambient silica dryer (design choice): the adsorption heat of the
    # silica gel is rejected to ambient air in a dry cooler, the feed enters the contactor at T0 + 5 K
    rej_ac = np.zeros_like(w['T0'])
    if case.drying == 'cond_silica' and case.cold == 'ambient':
        T_hot = L['T_ads']
        T_new = np.minimum(T_hot, w['T0'] + 5.0)
        rej_ac = plant.CP_AIR * (T_hot - T_new)                             # kW per kg/s air
        aktif = rej_ac > 0
        fan_ac = np.where(aktif, p.dP_coil / (P / (287.05 * w['T0'])) / p.eta_fan / 1e3, 0.0)
        el_ac = p.fan_rej_frac * rej_ac + fan_ac
        L['T_ads'] = T_new
        L['elec'] = L['elec'] + el_ac
        L['el_cold'] = L['el_cold'] + p.fan_rej_frac * rej_ac
        L['fan'] = L['fan'] + fan_ac
    L['rej_ac'] = rej_ac
    T_ads = L['T_ads']
    RH = np.clip(L['w_feed'] / humidity_ratio(T_ads, P), 0.0, 1.0)
    dq = pt('dq_CO2', T_ads, RH, log=True)                    # mol/kg
    purity = np.clip(pt('murni_pct', T_ads, RH), 0, 100) / 100.0
    n_rel = np.maximum(pt('n_rel', T_ads, RH), 1.0)
    w_rel = np.maximum(pt('w_rel', T_ads, RH), 0.0)
    qdes = pt('qdes', T_ads, RH)                               # kJ/mol CO2
    wcap = dq * M_CO2 * 1e-3                                   # kg CO2 / kg sorbent
    cap_ok = wcap > 0.1e-3 * M_CO2
    wcap_s = np.where(cap_ok, wcap, 1.0)
    x_co2 = p.x_co2_ppm * 1e-6 * M_CO2 / 28.96
    co2 = p.eta_cap * x_co2 * np.ones_like(T_ads)              # purity reported, not enforced
    c_bed = sb.cp + p.cp_contactor
    q_des = co2 * (c_bed / wcap_s * np.maximum(sb.T_des - np.maximum(w['T0'], T_ads), 0.0)
                   + qdes / (M_CO2 * 1e-3))
    V = tvsa.void_volume(sb, p)                                # m3 / kg sorbent
    P0, P1 = 101325.0, sb.p_des * 1e5
    w_evak = V * (P0 - P1 - P1 * np.log(P0 / P1)) / 1e3 / p.eta_vac        # kJ / kg sorbent / cycle
    vac = (co2 / M_CO2 * n_rel * 8.314 * sb.T_des * np.log(1.0 / sb.p_des) / p.eta_vac
           + co2 / wcap_s * w_evak)
    m_w = co2 * w_rel * M_H2O / M_CO2
    L.update(
        wcap=wcap, purity=purity, co2=co2,
        heat=L['heat'] - L['q_des'] + q_des, q_des=q_des,
        elec=L['elec'] - L['vac'] + vac, vac=vac,
        water=L['water'] + m_w, m_w=m_w, r_w=np.zeros_like(T_ads),
        cycles=plant.cycles_per_h(sb, p, wcap_s, np.zeros_like(T_ads), T_ads, c_bed),
        feasible=L['feasible'] & cap_ok,
    )
    return L


plant.air_train = air_train_gcmc
_CAPEX_ASLI = plant.process_capex


def process_capex_gcmc(L, op, case, sb):
    """dacsys capital cost plus the dry cooler of the after-cooled silica dryer (c_rej per kW rejected) and the
    replacement of the sorbent every UMUR_SORBEN years (present value; contactor structure not replaced; the
    replacements are scaled so that fixed O&M is charged on the initial capital only)."""
    p = case.prm
    parts = _CAPEX_ASLI(L, op, case, sb)
    extra = p.c_rej * plant._peak(L.get('rej_ac', 0.0 * L['co2']), op)
    parts['cold'] += extra
    parts['total'] += extra
    awal = sb.cost * p.sorbent_cost_factor * p.sorbent_fab * parts['sorbent_kg']
    ganti = awal * (faktor_ganti(p) - 1.0) * p.crf / (p.crf + p.om_proc)
    parts['sorbent'] += ganti
    parts['total'] += ganti
    return parts


plant.process_capex = process_capex_gcmc
PRM = Prm(purity_model='henry', p_des=0.01)


def kasus(site, koreksi=False):
    out = []
    for name, s in daftar_sorben(koreksi).items():
        drying = 'none' if s.hydrophobic else 'cond_silica'
        out.append(Case(site=site, drying=drying, cold='ambient', sorbent=name, lever='GCMC', prm=PRM))
    return out


def jalankan(site):
    w = cs.load(site)
    res = plant.resources(w, PRM)
    rows = []
    for kor in (False, True):
        for c in kasus(site, kor):
            if kor and ':kor' not in c.sorbent:
                continue
            r = plant.evaluate(c, w, res)
            rows.append(dict(site=site, sorbent=c.sorbent, drying=c.drying, **{k: v for k, v in r.items() if np.isscalar(v)}))
    return rows


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'uji'
    if mode == 'uji':
        import time
        t = time.time()
        rows = jalankan('Reykjavik')
        pd.set_option('display.width', 250)
        cols = ['site', 'sorbent', 'status', 'lcoc', 'hours', 'wcap_mmolg', 'purity', 'heat_GJ', 'heat_des_GJ', 'heat_sil_GJ',
                'elec_GJ', 'elec_vac_GJ', 'capex_sorbent', 'sorbent_kg', 'T_ads_mean']
        print(pd.DataFrame(rows).reindex(columns=cols).round(3).to_string(index=False))
        print('detik', round(time.time() - t, 1))
    else:
        from multiprocessing import Pool
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        with Pool(n) as pool:
            rows = sum(pool.map(jalankan, list(cs.SITES)), [])
        out = pd.DataFrame(rows)
        iklim = pd.read_csv(os.path.join(D, 'lcoc_amina.csv'))[['site', 'T_mean_C', 'RH_mean']]   # site climate statistics (lcoc_amina.py)
        out = out.merge(iklim, on='site', how='left')
        out.to_csv(os.path.join(D, 'lcoc_gcmc.csv' if PERLAKUAN == 'lima' else f'lcoc_gcmc_{PERLAKUAN}.csv'), index=False)
        print(out[['site', 'sorbent', 'status', 'lcoc', 'wcap_mmolg', 'purity', 'hours']].round(1).to_string(index=False))
