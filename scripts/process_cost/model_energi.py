"""Level-1 process model: equilibrium energy of one cycle and a lower-bound cost per tonne of CO2.

No dynamic column model: no mass/heat transfer, cycle time or equipment capital cost.
Scenario A : hydrophobic sorbent (KEN family), air at 50 % RH at T_ads, no drying.
Scenario B : strongly adsorbing hydrophilic sorbent; the air is first dried to the RH tolerated by the material.
Sub-ambient = cold-climate air as found (no refrigeration); the bed is cooled by the air itself.

usage: python model_energi.py            -> data/energi_ekonomi.csv, data/energi_terbaik.csv
input: data/kerja_final.csv (iterasi_kerja.py), kerja_hasil.csv, entalpi_ads298.csv,
       gcmcriil_final.csv, cp_semua.csv, saring_tahap1.csv
"""
import os
import numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
b = lambda f: pd.read_csv(os.path.join(D, f))
G = ['N2', 'O2', 'AR', 'H2O', 'CO2']
R, P0, M_CO2 = 8.314, 101325.0, 44.01          # J/mol/K, Pa, g/mol
KJMOL_KE_GJT = 1000 / M_CO2 / 1000             # 1 kJ/mol CO2 = 0.02272 GJ/t
PSAT = {298: 3170.0, 273: 611.0, 253: 103.0}   # Pa, same values as in the GCMC runs
Y_CO2 = 0.000415

# ---- parameters: (base, low, high) -------------------------------------------------
PAR = dict(
    harga_panas=(3.6, 3.0, 10.0),      # $/GJ; EIA industrial natural gas (2024) .. expensive heat
    harga_listrik=(0.087, 0.05, 0.15), # $/kWh; EIA industrial 2026
    ongkos_produksi=(50.0, 13.0, 71.0),# $/kg MOF excluding metal; DeSantis et al. 2017
    umur_siklus=(1e4, 1e3, 1e5),       # cycles; no data -> scanned
    pulih_panas=(0.0, 0.0, 0.5),       # fraction of sensible heat recovered
    eta_vakum=(0.5, 0.3, 0.7),         # isothermal vacuum pump efficiency
    dP_kipas=(200.0, 100.0, 500.0),    # Pa, contactor pressure drop (assumed)
    eta_kipas=(0.6, 0.6, 0.6),
    fraksi_tangkap=(0.75, 0.5, 0.9),   # fraction of the CO2 in the air that is captured
    panas_kering=(55.0, 45.0, 100.0),  # kJ/mol water removed (desiccant, thermal)
    rongga_unggun=(0.4, 0.4, 0.4),
)
HARGA_LOGAM = {'Ag': 2060.0,   # $/kg; ~64 $/oz, Trading Economics/Forbes 23-24 Sep 2026
               'Ni': 16.5,     # LME 23 Sep 2026
               'Nb': 57.0 / 0.699,  # from Nb2O5 at 57 $/kg (SMM Aug 2026), Nb mass fraction 0.699
               'Cd': 3.9,      # USGS MCS 2026
               'Zn': 4.0}      # LME Sep 2026 (metalcharts)
KOMPOSISI = {  # atoms per unit cell from the CIF
    'KEN':   {'Ag': 6, 'Br': 12, 'C': 144, 'Cl': 6, 'H': 108, 'N': 40, 'O': 24},
    'KEPFIX02_ion': {'Ag': 6, 'C': 144, 'Cl': 18, 'H': 108, 'N': 40, 'O': 24},
    'OCEGOV_ASR': {'C': 8, 'F': 6, 'H': 8, 'N': 4, 'Ni': 1, 'Si': 1},
    'ARAHIM01_ASR': {'C': 16, 'F': 12, 'H': 16, 'N': 8, 'Nb': 2, 'Ni': 2},
    'FUDTOQ_ASR': {'Br': 4, 'C': 76, 'Cd': 4, 'H': 52, 'N': 8, 'O': 16},
    'jp9b03221_si_002_ASR': {'C': 32, 'F': 24, 'H': 32, 'N': 16, 'Si': 4, 'Zn': 4},
}
MASSA = {'Ag': 107.868, 'Br': 79.904, 'C': 12.011, 'Cl': 35.453, 'H': 1.008, 'N': 14.007,
         'O': 15.999, 'F': 18.998, 'Ni': 58.693, 'Si': 28.086, 'Nb': 92.906, 'Cd': 112.41, 'Zn': 65.38}


def komposisi(m):
    return KOMPOSISI.get(m, KOMPOSISI['KEN'] if m.startswith('KEN') else None)


def ongkos_logam(m):
    k = komposisi(m); tot = sum(MASSA[e] * n for e, n in k.items())
    return sum(MASSA[e] * n / tot * HARGA_LOGAM[e] for e, n in k.items() if e in HARGA_LOGAM)


# ---- material data -------------------------------------------------------------------------
cp = b('cp_semua.csv').set_index('nama')
geo = b('saring_tahap1.csv'); geo['nama'] = geo.nama.str.replace('_pacman', '', regex=False)
geo = geo.set_index('nama')
hasil = b('kerja_hasil.csv')
fin = b('gcmcriil_final.csv').set_index(['nama', 'rh'])
h298 = b('entalpi_ads298.csv').pivot_table(index='nama', columns='komp', values='dH_kJmol')
hsub = hasil[hasil.run.str.contains('__ads')].assign(
    nama=lambda x: x.run.str.split('__').str[0], k=lambda x: x.run.str.split('__').str[1])


def entalpi(m, Ta):
    """|dH| of adsorption per component (kJ/mol) in the state after adsorption."""
    if Ta == 298:
        v = h298.loc[m]
    else:
        rh = 'rh050' if m.endswith('_ion') else 'rh000'
        v = hsub[(hsub.nama == m) & (hsub.k == f'ads{Ta}_{rh}')].set_index('komp').dH_kJmol
    return v.reindex(G).abs().fillna(0).values


def panas_sensibel(m, Ta, Td):
    """kJ/kg, integral of a linear Cp fitted to Cp at 300/350/400 K."""
    c = cp.loc[m]; a1, a0 = np.polyfit([300, 350, 400], [c.Cp_300K, c.Cp_350K, c.Cp_400K], 1)
    return a0 * (Td - Ta) + a1 / 2 * (Td ** 2 - Ta ** 2)       # J/g/K * K = kJ/kg


def geometri(m):
    g = geo.loc[m] if m in geo.index else geo[geo.basis == m.split('_')[0][:6]].iloc[0]
    return float(g.PV), float(g.rho)                            # cm3/g, g/cm3


# humidity tolerance in scenario B: largest RH that still retains >= 80 % of the dry uptake
rh_rendah = hasil[hasil.run.str.contains(r'__ads\d+_rh0(0[125]|10)') & (hasil.komp == 'CO2')].assign(
    nama=lambda x: x.run.str.split('__').str[0],
    Ta=lambda x: x.run.str.extract(r'ads(\d+)_')[0].astype(int),
    RH=lambda x: x.run.str.extract(r'_rh(\d+)$')[0].astype(int))


def rh_toleransi(m, Ta):
    kering = fin.loc[(m, 'rh000')].CO2 if Ta == 298 else \
        hsub[(hsub.nama == m) & (hsub.k == f'ads{Ta}_rh000') & (hsub.komp == 'CO2')].q_molkg.iat[0]
    x = rh_rendah[(rh_rendah.nama == m) & (rh_rendah.Ta == Ta)].sort_values('RH')
    ok = 0
    for _, r in x.iterrows():
        if r.q_molkg >= 0.8 * kering: ok = r.RH
        else: break
    return ok


def baris_toleransi(x):
    """Scenario B with a tolerated humidity > 0: the working capacity refers to adsorption from the feed dried to that
    humidity (GCMC adsorption state at RH_tol, desorption at 1 kPa from the surrogate), so that capacity and drying
    energy refer to the same feed. Returns the modified row, the row unchanged if RH_tol = 0, or None if no surrogate
    is available (P_des other than 1 kPa)."""
    m, Ta, Td = x.nama, int(x.Ta), int(x.Td)
    if m.endswith('_ion'):
        return x
    rh = rh_toleransi(m, Ta)
    if rh == 0:
        return x
    if x.P != 'p01' or Td not in (348, 373, 398):
        return None
    import surrogat_desorpsi as sd
    run = f'{m}__ads{Ta}_rh{rh:03d}'
    q_ads = hasil[hasil.run == run].set_index('komp').q_molkg.reindex(G).fillna(0).values
    dq, _ = sd.swakonsisten(q_ads, sd.k_dari(sd.KT, m, Td))
    y = x.copy()
    for g, v in zip(G, dq):
        y['d' + g] = v
    return y


def hitung(x, p):
    m, Ta, Td, Pd = x.nama, int(x.Ta), int(x.Td), {'p01': 1000.0, 'p05': 5000.0, 'p20': 20265.0}[x.P]
    dq = np.array([x['d' + g] for g in G]); dco2 = dq[4]
    if dco2 <= 0:
        return None
    lepas = np.clip(dq, 0, None)
    pv, rho = geometri(m)
    # heat (kJ/kg sorbent)
    q_sens = panas_sensibel(m, Ta, Td) * (1 - p['pulih_panas'])
    q_des = float(lepas @ entalpi(m, Ta))
    # vacuum electricity (kJ/kg): evacuation of the void gas, then compression of the product to 1 atm
    V = (pv + p['rongga_unggun'] / (1 - p['rongga_unggun']) / rho) * 1e-3     # m3/kg
    w_evak = V * (P0 - Pd - Pd * np.log(P0 / Pd)) / 1000 / p['eta_vakum']
    w_prod = lepas.sum() * R * 298 * np.log(P0 / Pd) / 1000 / p['eta_vakum']
    # fan and drying (per mol CO2)
    n_udara = 1 / (Y_CO2 * p['fraksi_tangkap'])
    w_kipas = n_udara * R * Ta / P0 * p['dP_kipas'] / p['eta_kipas'] / 1000       # kJ/mol CO2
    skenario = 'A' if m.endswith('_ion') else 'B'
    rh_tol = 50 if skenario == 'A' else rh_toleransi(m, Ta)
    n_air = n_udara * max(0.0, (50 - rh_tol) / 100 * PSAT[Ta] / P0)
    e_th = (q_sens + q_des) / dco2 + n_air * p['panas_kering']                     # kJ/mol CO2
    e_el = (w_evak + w_prod) / dco2 + w_kipas
    # cost ($/t CO2)
    gj = lambda e: e * KJMOL_KE_GJT
    c_energi = gj(e_th) * p['harga_panas'] + gj(e_el) * 277.8 * p['harga_listrik']
    c_sorben_kg = ongkos_logam(m) + p['ongkos_produksi']
    c_sorben = c_sorben_kg / (dco2 * M_CO2 / 1000 * p['umur_siklus']) * 1000
    murni = lepas[4] / (lepas.sum() - lepas[3]) * 100          # CO2 purity on a dry basis (water condensed)
    return dict(skenario=skenario, rh_toleransi=rh_tol, dq_CO2=dco2, murni_pct=murni, murni_basah_pct=lepas[4] / lepas.sum() * 100,
                panas_sensibel=gj(q_sens / dco2), panas_desorpsi=gj(q_des / dco2),
                pengeringan=gj(n_air * p['panas_kering']), listrik=gj(e_el),
                energi_GJt=gj(e_th + e_el), sorben_kg=c_sorben_kg,
                ongkos_energi=c_energi, ongkos_sorben=c_sorben, ongkos_total=c_energi + c_sorben)


if __name__ == '__main__':
    kf = b('kerja_final.csv')
    dasar = {k: v[0] for k, v in PAR.items()}
    out = []
    for _, x in kf.iterrows():
        x = baris_toleransi(x)
        if x is None:
            continue
        h = hitung(x, dasar)
        if h: out.append({**x[['nama', 'Ta', 'Td', 'P', 'it', 'stabil']].to_dict(), **h})
    e = pd.DataFrame(out)
    e.to_csv(os.path.join(D, 'energi_ekonomi.csv'), index=False)
    # cheapest condition per material and T_ads, then one-at-a-time sensitivity
    best = e.sort_values('ongkos_total').groupby(['nama', 'Ta']).head(1).sort_values(['Ta', 'ongkos_total'])
    sens = []
    for _, x0 in best.iterrows():
        x = baris_toleransi(kf[(kf.nama == x0.nama) & (kf.Ta == x0.Ta) & (kf.Td == x0.Td) & (kf.P == x0.P)].iloc[0])
        for k, (v0, lo, hi) in PAR.items():
            if lo == hi: continue
            r = [hitung(x, {**dasar, k: v})['ongkos_total'] for v in (lo, hi)]
            sens.append(dict(nama=x0.nama, Ta=x0.Ta, parameter=k, rendah=r[0], tinggi=r[1], dasar=x0.ongkos_total))
    best.to_csv(os.path.join(D, 'energi_terbaik.csv'), index=False)
    pd.DataFrame(sens).to_csv(os.path.join(D, 'energi_sensitivitas.csv'), index=False)
    pd.set_option('display.width', 250)
    kol = ['nama', 'Ta', 'Td', 'P', 'stabil', 'skenario', 'rh_toleransi', 'dq_CO2', 'murni_pct', 'panas_sensibel',
           'panas_desorpsi', 'pengeringan', 'listrik', 'energi_GJt', 'sorben_kg', 'ongkos_energi', 'ongkos_sorben', 'ongkos_total']
    print(best[kol].round(2).to_string(index=False))
