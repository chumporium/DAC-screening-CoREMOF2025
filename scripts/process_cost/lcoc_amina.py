"""Amine DAC reference on the same off-grid energy block and sites as dacsys_gcmc.py.

The dacsys package (path in the environment variable DACSYS_PATH) is used read-only (no bytecode, output to this project).
Energy demand of amine DAC from the fit to Wenzel et al. (dacsys/amine_benchmark.json) at the annual mean
temperature and humidity of each site; capital cost 600/1000/1500 USD per (t/yr), O&M 3/5/10 %, energy
factor 0.75/1.0/1.25 (lo/base/hi), all dacsys defaults.

python lcoc_amina.py [n_proc] -> data/lcoc_amina.csv
"""
import os, sys
sys.dont_write_bytecode = True
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
ARS = os.environ.get('DACSYS_PATH', 'dacsys_path_not_set')      # folder that contains the dacsys package
sys.path.insert(0, ARS)

import json
import pandas as pd
from dacsys import plant, climate_sites as cs
from dacsys.levers import Prm

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(ROOT, 'data')


def jalankan(site):
    w = cs.load(site)
    r = plant.amine(w, p=Prm())
    return dict(site=site, **r)


if __name__ == '__main__':
    from multiprocessing import Pool
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    with Pool(n) as pool:
        rows = pool.map(jalankan, list(cs.SITES))
    out = pd.DataFrame(rows)
    b = json.load(open(os.path.join(ARS, 'dacsys', 'amine_benchmark.json')))
    out['within_fit_range'] = (out.T_mean_C.between(*b['T_range']) & out.RH_mean.between(*b['RH_range']))
    out.to_csv(os.path.join(D, 'lcoc_amina.csv'), index=False)
    pd.set_option('display.width', 200)
    print(out.round(1).to_string(index=False))
