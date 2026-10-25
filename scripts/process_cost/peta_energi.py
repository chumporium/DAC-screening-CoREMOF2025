"""Working capacity, purity and equilibrium energy on the adsorption grid (T_ads x RH), P_des = 1 kPa.

Adsorption states: five-component GCMC (grid + earlier runs). Desorption states: surrogate
calibrated on the self-consistent GCMC runs (surrogat_desorpsi.py). Energy terms as in model_energi.py.
python peta_energi.py -> data/peta_iklim.csv
Simplified treatments for the system-level comparison (environment variable PERLAKUAN):
  inert  N2, O2 and Ar neither adsorbed nor released (self-consistent CO2/H2O release) -> peta_iklim_inert.csv
  umpan  preset desorption gas with the composition of the feed at P_des            -> peta_iklim_umpan.csv
"""
import os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_energi as me
import surrogat_desorpsi as sd

G = me.G
PERLAKUAN = os.environ.get('PERLAKUAN', 'lima')
Y_UDARA = np.array([0.78084, 0.20946, 0.00934, 0.0, 4.15e-4])   # N2, O2, Ar, H2O, CO2 (dry air)
dasar = {k: v[0] for k, v in me.PAR.items()}
hasil = me.b('kerja_hasil.csv')
hasil = hasil[hasil.run.str.match(r'.+__ads\d+_rh\d+$')]
ads = hasil.pivot_table(index='run', columns='komp', values='q_molkg').reindex(columns=G).fillna(0)
dH = hasil.pivot_table(index='run', columns='komp', values='dH_kJmol').reindex(columns=G)
fin = me.b('gcmcriil_final.csv')
h298 = me.b('entalpi_ads298.csv').pivot_table(index='nama', columns='komp', values='dH_kJmol').reindex(columns=G)
PSAT = {298: 3170.0, 273: 611.0, 258: 163.0, 253: 102.0, 243: 37.4}

titik = []   # (nama, Ta, RH, q_ads[G], |dH|[G])
for run in ads.index:
    m = re.match(r'(.+)__ads(\d+)_rh(\d+)$', run)
    nama, Ta, RH = m.group(1), int(m.group(2)), int(m.group(3))
    titik.append((nama, Ta, RH, ads.loc[run].values, dH.loc[run].abs().values))
for _, r in fin.iterrows():   # 298 K states from the screening GCMC
    if r.nama in me.cp.index and r.nama in sd.KT.nama.values:
        RH = int(r.rh[2:])
        h = (h298.loc[r.nama].abs().values if (r.nama in h298.index and ((RH == 50) == r.nama.endswith('_ion') or RH == 0 and not r.nama.endswith('_ion')))
             else np.full(len(G), np.nan))
        titik.append((r.nama, 298, RH, r[G].astype(float).fillna(0).values, h))

rows = []
for nama, Ta, RH, q_ads, h in titik:
    if nama not in sd.KT.nama.values or q_ads[4] <= 0:
        continue
    if nama.endswith('_ion') is False and RH > 10:
        continue                                     # scenario B only with a dried feed
    h = np.where(np.isfinite(h), h, np.nanmedian(h) if np.isfinite(h).any() else 0.0)
    pv, rho = me.geometri(nama)
    for Td in (348, 373, 398):
        k = sd.k_dari(sd.KT, nama, Td)
        if not k.any():
            continue
        if PERLAKUAN == 'inert':
            qa = q_ads.copy(); qa[:3] = 0.0
            dq, y = sd.swakonsisten(qa, k)
        elif PERLAKUAN == 'umpan':
            yh = PSAT[Ta] * RH / 100.0 / me.P0 if Ta in PSAT else 0.0
            yf = Y_UDARA * (1.0 - yh); yf[3] = yh
            dq = q_ads - k * yf * sd.P_DES; y = yf
        else:
            dq, y = sd.swakonsisten(q_ads, k)
        if dq[4] <= 0:
            rows.append(dict(nama=nama, Ta=Ta, RH=RH, Td=Td, dq_CO2=dq[4], murni_pct=np.nan)); continue
        lepas = np.clip(dq, 0, None)
        q_sens = me.panas_sensibel(nama, Ta, Td)
        q_des = float(lepas @ h)
        V = (pv + dasar['rongga_unggun'] / (1 - dasar['rongga_unggun']) / rho) * 1e-3
        Pd = sd.P_DES
        w_el = (V * (me.P0 - Pd - Pd * np.log(me.P0 / Pd)) / 1000 + lepas.sum() * me.R * 298 * np.log(me.P0 / Pd) / 1000) / dasar['eta_vakum']
        n_udara = 1 / (me.Y_CO2 * dasar['fraksi_tangkap'])
        w_kipas = n_udara * me.R * Ta / me.P0 * dasar['dP_kipas'] / dasar['eta_kipas'] / 1000
        gj = lambda e: e * me.KJMOL_KE_GJT
        rows.append(dict(nama=nama, Ta=Ta, RH=RH, Td=Td, dq_CO2=dq[4], murni_pct=lepas[4] / (lepas.sum() - lepas[3]) * 100,
                         murni_basah_pct=lepas[4] / lepas.sum() * 100,
                         H2O_lepas=lepas[3], N2O2Ar_lepas=lepas[0] + lepas[1] + lepas[2],
                         panas_sensibel=gj(q_sens / dq[4]), panas_desorpsi=gj(q_des / dq[4]),
                         listrik=gj(w_el / dq[4] + w_kipas), energi_tanpa_kering=gj((q_sens + q_des + w_el) / dq[4] + w_kipas)))
out = pd.DataFrame(rows).sort_values(['nama', 'Ta', 'RH', 'Td'])
out.to_csv(os.path.join(me.D, 'peta_iklim.csv' if PERLAKUAN == 'lima' else f'peta_iklim_{PERLAKUAN}.csv'), index=False)
pd.set_option('display.width', 220)
x = out[(out.Td == 398) & out.nama.isin(['KENVOR_ion', 'KENVOR02_ion'])]
print(x.pivot_table(index=['nama', 'Ta'], columns='RH', values='dq_CO2').round(3).to_string())
print(x.pivot_table(index=['nama', 'Ta'], columns='RH', values='murni_pct').round(0).to_string())
print(x.pivot_table(index=['nama', 'Ta'], columns='RH', values='energi_tanpa_kering').round(1).to_string())
y = out[(out.Td == 398) & ~out.nama.str.endswith('_ion') & (out.RH == 0)]
print(y.pivot_table(index='nama', columns='Ta', values='dq_CO2').round(3).to_string())
