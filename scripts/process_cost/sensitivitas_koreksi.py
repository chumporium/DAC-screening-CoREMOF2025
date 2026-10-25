"""Scenario-B working capacity and energy corrected with experiment/simulation uptake ratios.

f = q_exp / q_GCMC (298 K, dry, ~400-500 ppm); dq_CO2 is scaled by f, other components unchanged.
This assumes the residual loading after desorption is overestimated by the same factor.
python sensitivitas_koreksi.py -> data/sensitivitas_koreksi.csv
"""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_energi as me

F = {  # (q_exp mol/kg, condition, q_GCMC dry air 298 K)
    'ARAHIM01_ASR': (1.3, '400 ppm, Bhatt 2016', 1.78),
    'OCEGOV_ASR': (0.4, '500 ppm, Mukherjee 2019', 0.8157),
    'jp9b03221_si_002_ASR': (0.13, '400 ppm, Shekhah 2014 (SIFSIX-3-Zn)', 1.682),
    'FUDTOQ_ASR': (None, 'no experimental data', 0.9992),
}
kf = me.b('kerja_final.csv')
dasar = {k: v[0] for k, v in me.PAR.items()}
rows = []
for m, (qe, src, qs) in F.items():
    f = qe / qs if qe else 1.0
    for Ta in (298, 273, 253):
        x = me.baris_toleransi(kf[(kf.nama == m) & (kf.Ta == Ta) & (kf.Td == 398) & (kf.P == 'p01')].iloc[0]).copy()
        h0 = me.hitung(x, dasar)
        x['dCO2'] = x['dCO2'] * f
        h1 = me.hitung(x, dasar)
        rows.append(dict(nama=m, f=round(f, 3), sumber=src, Ta=Ta, dq=h0['dq_CO2'], dq_kor=h1['dq_CO2'],
                         E=h0['energi_GJt'], E_kor=h1['energi_GJt'], murni=h0['murni_pct'], murni_kor=h1['murni_pct']))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(me.D, 'sensitivitas_koreksi.csv'), index=False)
pd.set_option('display.width', 200)
print(out.round(2).to_string(index=False))
a = me.b('energi_ekonomi.csv')
print(a[(a.nama.str.endswith('_ion')) & (a.Td == 398) & (a.P == 'p01')].groupby('Ta')[['dq_CO2', 'energi_GJt']].agg(['min', 'max']).round(2))
