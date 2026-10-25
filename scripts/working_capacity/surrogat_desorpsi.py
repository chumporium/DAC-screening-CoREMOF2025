"""Desorption surrogate at P_des = 1 kPa: residual loading q_i = k_i(T_des) * y_i * P_des.

k_i is calibrated per material, component and T_des from all GCMC desorption runs at 1 kPa
(every iteration is a valid GCMC state). The desorption gas composition is then solved
self-consistently for any adsorption state, with the same damped iteration as the GCMC campaign.

python surrogat_desorpsi.py validasi         -> leave-one-T_ads-out test against GCMC (kerja_final), 1 kPa
python surrogat_desorpsi.py validasi p05     -> the same for 5 kPa (surrogat_validasi_p05.csv)
"""
import glob, os, re, sys
import numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
G = ['N2', 'O2', 'AR', 'H2O', 'CO2']
P_DES = 1000.0
TEKANAN = {'p01': 1000.0, 'p05': 5000.0}

hasil = pd.read_csv(os.path.join(D, 'kerja_hasil.csv'))
q = hasil.pivot_table(index='run', columns='komp', values='q_molkg').reindex(columns=G).fillna(0)
komp = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(D, 'compositions', 'komp_*.csv'))
                  if 'rh_rendah' not in f and 'komp_ads' not in f and 'grid' not in f]).drop_duplicates('run').set_index('run')


def kalibrasi_tabel(kode='p01'):
    P = TEKANAN[kode]; pola = re.compile(r'(.+)__a(\d+)_d(\d+)_' + kode + r'_i(\d+)$')
    rows = []
    for run in q.index:
        m = pola.match(run)
        if not m or run not in komp.index:
            continue
        nama, Ta, Td, it = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
        for i, g in enumerate(G):
            y = float(komp.loc[run, 'y' + g])
            if y > 1e-6 and q.loc[run, g] > 0:
                rows.append(dict(nama=nama, Ta=Ta, Td=Td, it=it, komp=g, k=q.loc[run, g] / (y * P)))
    return pd.DataFrame(rows)


KT = kalibrasi_tabel()


def k_dari(kt, nama, Td):
    x = kt[(kt.nama == nama) & (kt.Td == Td)].groupby('komp').k.median()
    return np.array([x.get(g, 0.0) for g in G])


def swakonsisten(q_ads, k, n_iter=60, tol=1e-4, P=P_DES):
    """q_ads: adsorbed amounts (mol/kg) in order G. Returns dq (mol/kg), y_des."""
    y = q_ads / q_ads.sum()
    for _ in range(n_iter):
        q_des = k * y * P
        dq = q_ads - q_des
        lepas = np.clip(dq, 0, None)
        if lepas.sum() <= 0:
            return dq, y
        y_baru = 0.5 * y + 0.5 * lepas / lepas.sum()
        if np.abs(y_baru - y).max() < tol:
            y = y_baru; break
        y = y_baru
    return q_ads - k * y * P, y


def validasi(kode='p01'):
    P = TEKANAN[kode]; KTk = KT if kode == 'p01' else kalibrasi_tabel(kode)
    kf = pd.read_csv(os.path.join(D, 'kerja_final.csv'))
    fin = pd.read_csv(os.path.join(D, 'gcmcriil_final.csv')).set_index(['nama', 'rh'])
    rows = []
    for _, x in kf[(kf.P == kode) & kf.Td.isin([348, 373, 398])].iterrows():
        m, Ta, Td = x.nama, int(x.Ta), int(x.Td)
        rh = 'rh050' if m.endswith('_ion') else 'rh000'
        if Ta == 298:
            q_ads = fin.loc[(m, rh)][G].astype(float).fillna(0).values
        else:
            q_ads = q.loc[f'{m}__ads{Ta}_{rh}'].values
        kt = KTk[KTk.Ta != Ta]                    # leave the adsorption temperature out
        dq, y = swakonsisten(q_ads, k_dari(kt, m, Td), P=P)
        rows.append(dict(nama=m, Ta=Ta, Td=Td, dq_gcmc=x.dq, sd=x.sd, dq_sur=dq[4], stabil=x.stabil))
    v = pd.DataFrame(rows)
    v['selisih'] = v.dq_sur - v.dq_gcmc
    v['rel_pct'] = v.selisih / v.dq_gcmc.abs().clip(lower=0.02) * 100
    v.to_csv(os.path.join(D, 'surrogat_validasi.csv' if kode == 'p01' else f'surrogat_validasi_{kode}.csv'), index=False)
    pd.set_option('display.width', 200)
    for Td in (348, 373, 398):
        s = v[v.Td == Td]
        print(f'T_des {Td}: n={len(s)}  median |selisih| {s.selisih.abs().median():.3f} mol/kg, '
              f'maks {s.selisih.abs().max():.3f}; median |rel| {s.rel_pct.abs().median():.1f} %, '
              f'dalam 2 sd: {(s.selisih.abs() <= 2 * s.sd).mean() * 100:.0f} %')
    print(v[v.Td == 398].round(3).to_string(index=False))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'validasi':
        validasi(sys.argv[2] if len(sys.argv) > 2 else 'p01')
