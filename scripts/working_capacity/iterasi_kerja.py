"""Self-consistent working capacity: compute dq for each desorption run, check whether the
iteration is stable, and write the composition table of the next iteration for conditions
that are not yet stable.

Stable    : |dq_CO2(n) - dq_CO2(n-1)| <= 2 * sqrt(sd_n^2 + sd_(n-1)^2)
            (the difference between two independent estimates is not statistically significant)
Iteration : y(n+1) = 0.5 y(n) + 0.5 * released(n) / sum(released(n)),  released = max(0, q_ads - q_des)

usage: python iterasi_kerja.py <out_composition.csv>
input (in data/): kerja_hasil.csv, gcmcriil_final.csv, compositions/komp_*.csv
output: data/kerja_dq.csv (all runs), data/kerja_final.csv (last iteration of each condition)
"""
import glob, os, re, sys
import numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
G = ['N2', 'O2', 'AR', 'H2O', 'CO2']
Y = ['y' + g for g in G]
POLA = re.compile(r'(.+)__a(\d+)_d(\d+)_(?:(p\d+)_)?i(\d+)$')

r = pd.read_csv(os.path.join(D, 'kerja_hasil.csv'))
q = r.pivot_table(index='run', columns='komp', values='q_molkg').reindex(columns=G).fillna(0)
sg = r[r.komp == 'CO2'].set_index('run').sigma_pct
fin = pd.read_csv(os.path.join(D, 'gcmcriil_final.csv')).set_index(['nama', 'rh'])
komp = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(D, 'compositions', 'komp_*.csv'))
                  if 'rh_rendah' not in f and 'ads' not in os.path.basename(f)]).drop_duplicates('run').set_index('run')


def rh_ads(m):
    return 'rh050' if m.endswith('_ion') else 'rh000'


def ads(m, Ta):
    """(q per component, absolute CO2 sd) of the state after adsorption."""
    if Ta == 298:
        x = fin.loc[(m, rh_ads(m))]
        return x[G].astype(float).fillna(0).values, x.CO2 * x.sig_CO2 / 100
    run = f'{m}__ads{Ta}_{rh_ads(m)}'
    return q.loc[run].values, q.loc[run, 'CO2'] * sg[run] / 100


rows = []
for run in q.index:
    k = POLA.match(run)
    if not k:
        continue
    m, Ta, Td, P, it = k.group(1), int(k.group(2)), int(k.group(3)), k.group(4) or 'p20', int(k.group(5))
    a, sa = ads(m, Ta)
    d = a - q.loc[run].values
    lepas = np.clip(d, 0, None)
    sdes = q.loc[run, 'CO2'] * sg[run] / 100
    rows.append(dict(run=run, nama=m, Ta=Ta, Td=Td, P=P, it=it, dq=d[4], sd=np.hypot(sa, sdes),
                     lepas_pct=d[4] / a[4] * 100, murni=lepas[4] / lepas.sum() * 100 if lepas.sum() > 0 else np.nan,
                     **{'d' + g: d[i] for i, g in enumerate(G)}))
dq = pd.DataFrame(rows)
dq.to_csv(os.path.join(D, 'kerja_dq.csv'), index=False)

# stability: compare the last iteration with the previous one
kunci = ['nama', 'Ta', 'Td', 'P']
akhir = dq.sort_values('it').groupby(kunci).tail(1).set_index(kunci)
seb = dq.sort_values('it').groupby(kunci).nth(-2).set_index(kunci)
akhir['beda'] = (akhir.dq - seb.dq.reindex(akhir.index)).abs()
akhir['batas'] = 2 * np.hypot(akhir.sd, seb.sd.reindex(akhir.index))
akhir['stabil'] = akhir.beda <= akhir.batas
akhir.reset_index().to_csv(os.path.join(D, 'kerja_final.csv'), index=False)

baru, tanpa = [], []
for (m, Ta, Td, P), x in akhir.iterrows():
    if x.stabil:
        continue
    run = x.run
    d = x[['d' + g for g in G]].astype(float).values
    lepas = np.clip(d, 0, None)
    if lepas.sum() <= 0:
        tanpa.append(run); continue
    y = 0.5 * komp.loc[run, Y].astype(float).values + 0.5 * lepas / lepas.sum()
    baru.append([re.sub(r'_i\d+$', f'_i{x.it + 1}', run), Td, komp.loc[run, 'P'], *y])

print(f'kondisi: {len(akhir)}  stabil: {int(akhir.stabil.sum())}  '
      f'belum dibandingkan (baru i1): {int(akhir.beda.isna().sum())}  iterasi baru: {len(baru)}')
if tanpa:
    print('tanpa pelepasan sama sekali (dihentikan):', tanpa)
if len(sys.argv) > 1:
    pd.DataFrame(baru, columns=['run', 'T', 'P'] + Y).to_csv(sys.argv[1], index=False)
    print('->', sys.argv[1])
