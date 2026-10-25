"""Working capacity at P_des = 1 kPa under four treatments of the desorption atmosphere (Section 3.6, Table 3).

  feed    desorption gas with the CO2 mole fraction of the feed (p_CO2 = 4.15e-4 x P_des): fixed state
  pureCO2 desorption gas of pure CO2 at P_des: fixed state
  binary  self-consistent, but only CO2 and H2O adsorb and are released (N2, O2, Ar inert)
  five    self-consistent with all five components (this work)

The residual loadings follow the surrogate q_i = k_i(T_des) y_i P_des (Section 2.7), calibrated on the
GCMC desorption runs at 1 kPa, so that all four treatments use the same molecular model.
Purity: CO2 mole fraction of the released gas on a dry basis (water removed).

python banding_desorpsi.py  ->  data/banding_desorpsi.csv
"""
import os
import numpy as np, pandas as pd
import surrogat_desorpsi as sd

D = sd.D; G = sd.G; P = 1000.0; Y_FEED = 4.15e-4
KEN = ['KENVOR_ion', 'KENVIL_ion', 'KENTEF05_ion', 'KENTEF_ion', 'KENVOR02_ion', 'KEPFIX02_ion']
BGRP = ['OCEGOV_ASR', 'ARAHIM01_ASR', 'FUDTOQ_ASR', 'jp9b03221_si_002_ASR']
fin = pd.read_csv(os.path.join(D, 'gcmcriil_final.csv')).set_index(['nama', 'rh'])
P20 = 20265.0
kd = pd.read_csv(os.path.join(D, 'kerja_dq.csv')).set_index(['nama', 'Ta', 'Td', 'P', 'it'])
kf = pd.read_csv(os.path.join(D, 'kerja_final.csv')).set_index(['nama', 'Ta', 'Td', 'P'])
lab = pd.read_csv(os.path.join(D, 'kerja_label_tanpa_kapasitas.csv')).set_index(['nama', 'Ta', 'Td', 'P'])


def rayleigh(q_ads, k, P=P, langkah=4000):
    """Differential (Rayleigh-type) desorption at constant P_des: the gas released at each instant is in equilibrium
    with the current loading (p_i = q_i / k_i) and is removed at once; desorption ends when sum_i q_i / k_i = P.
    Brackets the single-stage flash of Eq. (3), in which the whole released gas stays in contact with the bed."""
    q = q_ads.astype(float).copy(); aktif = k > 0
    dN = q[aktif].sum() / langkah
    for _ in range(100 * langkah):
        p = np.where(aktif, q / np.where(aktif, k, 1), 0.0)
        S = p.sum()
        if S <= P:
            break
        y = p / S
        q = np.maximum(q - y * min(dN, (S - P) / (p / np.maximum(q, 1e-30)).max()), 0.0)
    return q_ads - q


def murni_kering(dq):
    lepas = np.clip(dq, 0, None); tot = lepas.sum() - lepas[3]
    return 100 * lepas[4] / tot if tot > 0 else np.nan


rows = []
for m in KEN + BGRP:
    rh = 'rh050' if m in KEN else 'rh000'
    for Ta in (298, 253):
        q_ads = (fin.loc[(m, rh)][G].astype(float).fillna(0).values if Ta == 298
                 else sd.q.loc[f'{m}__ads{Ta}_{rh}'].values)
        for Td in (373, 398):
            k = sd.k_dari(sd.KT, m, Td)
            dq5, y5 = sd.swakonsisten(q_ads, k)
            qb = q_ads.copy(); qb[:3] = 0.0
            dq2, y2 = sd.swakonsisten(qb, k)
            dqr = rayleigh(q_ads, k)
            rows.append(dict(nama=m, skenario='A' if m in KEN else 'B', Ta=Ta, Td=Td, q_ads_CO2=q_ads[4],
                             dq_feed=q_ads[4] - k[4] * Y_FEED * P, dq_pureCO2=q_ads[4] - k[4] * P,
                             dq_binary=dq2[4], dq_five=dq5[4], dq_rayleigh=dqr[4], purity_rayleigh=murni_kering(dqr),
                             purity_binary=murni_kering(dq2), purity_five=murni_kering(dq5),
                             yCO2_des_five=y5[4], yCO2_des_binary=y2[4],
                             # 20.3 kPa: feed-composition state (Henry regime, p_CO2 = 8.4 Pa) and GCMC
                             dq_feed_p20=q_ads[4] - k[4] * Y_FEED * P20,
                             dq_first_p20=kd.dq.get((m, Ta, Td, 'p20', 1), np.nan),
                             dq_five_p20=kf.dq.get((m, Ta, Td, 'p20'), np.nan),
                             sd_five_p20=kf.sd.get((m, Ta, Td, 'p20'), np.nan),
                             label_p20=lab.label.get((m, Ta, Td, 'p20'), '')))
r = pd.DataFrame(rows)
for c in ('feed', 'pureCO2', 'binary', 'rayleigh'):
    r[f'ratio_{c}'] = r[f'dq_{c}'] / r.dq_five
r.to_csv(os.path.join(D, 'banding_desorpsi.csv'), index=False)
pd.set_option('display.width', 220)
print(r.round(3).to_string(index=False))
for s in ('A', 'B'):
    for Ta in (298, 253):
        x = r[(r.skenario == s) & (r.Ta == Ta)]
        print(s, Ta, {c: (round(x[f'ratio_{c}'].min(), 2), round(x[f'ratio_{c}'].max(), 2)) for c in ('feed', 'pureCO2', 'binary', 'rayleigh')},
              'purity five', round(x.purity_five.min()), '-', round(x.purity_five.max()))
print(r[['nama', 'Ta', 'Td', 'dq_feed_p20', 'dq_first_p20', 'dq_five_p20', 'sd_five_p20', 'label_p20']].round(3).to_string(index=False))
