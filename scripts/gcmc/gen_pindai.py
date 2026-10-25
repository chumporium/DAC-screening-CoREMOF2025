"""RASPA3 input generator for five-component GCMC runs (N2, O2, Ar, H2O, CO2) at a given
relative humidity and, optionally, a modified CO2 mole fraction.

Framework and adsorbate models as in gen_gcmc.py (UFF + PACMAN-DDEC6 charges; TraPPE, TIP4P).
tulis() writes one run directory (CIF, force field, molecule and simulation JSON files) and is
also used by bangun_kerja.py. With grid=True, per-component density grids (.cube) and
conventional RDFs are written and all components use standard swap moves: RASPA3 samples
density grids only for whole molecules, so CFCMC (fractional CO2 molecules) would leave the
CO2 grid empty. With grid=False, CO2 uses CFCMC swap moves for better statistics.
RH 50 % corresponds to p_H2O = 1585 Pa (p_sat = 3170 Pa).

usage: python3 gen_pindai.py <dir_cif> <dir_out> <RH_percent|semua>
                             [cycles] [T_K] [x_CO2|-] [grid|nogrid]
"""
import sys, os, json, math, csv, io

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_widom import UFF, MASSA, baca_cif, lebar_tegak
from gen_gcmc import GAS, X_KERING, P_TOTAL

# water saturation pressure used for all humid runs: RH 50 % -> p_H2O = 1585 Pa
P_SAT = 3170.0
RH_DAFTAR = [0, 10, 25, 50, 75, 90]
KISI_RAPAT = [50, 50, 50]      # .cube grid per component; 50^3 points keeps files small


def fraksi_rh(rh, x_co2=None):
    """Mole fractions at relative humidity rh (percent).

    x_co2 overrides the CO2 mole fraction of air (0.000415); the other dry
    gases are rescaled so that the total is unchanged. Used for diagnostic runs only.
    """
    ker = dict(X_KERING)
    if x_co2 is not None:
        lain = {k: v for k, v in ker.items() if k != 'CO2'}
        skala = (1.0 - x_co2) / sum(lain.values())
        ker = {k: v * skala for k, v in lain.items()}
        ker['CO2'] = x_co2
    if rh <= 0:
        return ker
    x_air = rh / 100.0 * P_SAT / P_TOTAL
    sisa = 1.0 - x_air
    d = {k: v * sisa for k, v in ker.items()}
    d['H2O'] = x_air
    return d


def tulis(cif, keluar, rh, siklus=50000, suhu=298.0, cutoff=9.0,
          x_co2=None, grid=True):
    nama = os.path.basename(cif)[:-4]
    sel, atom = baca_cif(cif)
    if not atom:
        raise ValueError('CIF tanpa atom: %s' % cif)
    hilang = sorted({a['el'] for a in atom} - set(UFF))
    if hilang:
        raise ValueError('tanpa UFF: ' + ','.join(hilang))
    os.makedirs(keluar, exist_ok=True)

    tipe = {}
    for a in atom:
        k = (a['el'], round(float(a['q']), 3))
        if k not in tipe:
            tipe[k] = 'T%d' % len(tipe)
    q = "'"
    with open(os.path.join(keluar, nama + '.cif'), 'w', newline='\n') as f:
        f.write('data_%s\n' % nama)
        for tag, v in zip(['_cell_length_a', '_cell_length_b', '_cell_length_c',
                           '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'], sel):
            f.write('%-22s %s\n' % (tag, v))
        f.write('\n_symmetry_space_group_name_H-M    "P 1"\n')
        f.write('_symmetry_Int_Tables_number       1\n\n')
        f.write('loop_\n_symmetry_equiv_pos_as_xyz\n  %sx, y, z%s\n\n' % (q, q))
        f.write('loop_\n  _atom_site_type_symbol\n  _atom_site_label\n')
        f.write('  _atom_site_symmetry_multiplicity\n')
        f.write('  _atom_site_fract_x\n  _atom_site_fract_y\n  _atom_site_fract_z\n')
        f.write('  _atom_site_occupancy\n')
        for a in atom:
            f.write('%-3s %-8s 1.0  %s %s %s %s\n' % (
                a['el'], tipe[(a['el'], round(float(a['q']), 3))],
                a['x'], a['y'], a['z'], a['occ']))

    ps, si = [], []
    for (el, q3), nm in sorted(tipe.items(), key=lambda z: int(z[1][1:])):
        ps.append({"name": nm, "framework": True, "print_to_output": False, "element": el,
                   "print_as": el, "mass": MASSA.get(el, 1.0), "charge": q3,
                   "source": "PACMAN-DDEC6"})
        si.append({"name": nm, "type": "lennard-jones", "parameters": list(UFF[el]),
                   "source": "UFF"})
    x = fraksi_rh(rh, x_co2)
    sudah = set()
    for gas in x:
        m = GAS[gas]
        for (nm, sg, ep, qq, z, ms) in m['situs']:
            if nm in sudah:
                continue
            sudah.add(nm)
            ps.append({"name": nm, "framework": False, "print_to_output": sg is not None,
                       "element": nm[0], "print_as": nm[0] if sg is not None else '-',
                       "mass": ms, "charge": qq, "source": m['sumber']})
            si.append({"name": nm, "type": "lennard-jones" if sg is not None else "none",
                       "parameters": [ep, sg] if sg is not None else [0.0, 1.0],
                       "source": m['sumber']})
    json.dump({"MixingRule": "Lorentz-Berthelot", "TruncationMethod": "shifted",
               "TailCorrections": False, "PseudoAtoms": ps, "SelfInteractions": si},
              open(os.path.join(keluar, 'force_field.json'), 'w'), indent=1)

    for gas in x:
        m = GAS[gas]
        situs = ([[s[0], k] for s, k in zip(m['situs'], m['koord'])] if 'koord' in m
                 else [[s[0], [0.0, 0.0, s[4]]] for s in m['situs']])
        json.dump({"CriticalTemperature": m['Tc'], "CriticalPressure": m['Pc'],
                   "AcentricFactor": m['w'], "Type": "rigid",
                   "pseudoAtoms": situs, "Bonds": m['ikatan']},
                  open(os.path.join(keluar, gas + '.json'), 'w'), indent=1)

    komp = []
    for gas, xf in sorted(x.items(), key=lambda z: -z[1]):
        c = {"Name": gas, "MolFraction": round(xf, 10), "FugacityCoefficient": 1.0,
             "TranslationProbability": 0.5, "ReinsertionProbability": 0.5,
             "WidomProbability": 0.2, "CreateNumberOfMolecules": 0}
        if gas != 'AR':
            c["RotationProbability"] = 0.5
        # grid on -> standard swap moves for all components, because CFCMC prevents
        # density-grid sampling. grid off -> CFCMC for CO2, better statistics.
        if gas == 'CO2' and not grid:
            c["CFCMC_CBMC_SwapProbability"] = 1.0
            c["ThermodynamicIntegration"] = True
        else:
            c["SwapProbability"] = 1.0
        komp.append(c)

    n = [max(1, int(math.ceil(2 * cutoff / w))) for w in lebar_tegak(sel)]
    sistem = {"Type": "Framework", "Name": nama, "NumberOfUnitCells": n,
              "ExternalTemperature": suhu, "ExternalPressure": P_TOTAL,
              "ChargeMethod": "Ewald", "CutOffVDW": cutoff, "CutOffCoulomb": cutoff,
              # density grids and RDFs show where each adsorbate is located
              "ComputeDensityGrid": bool(grid),
              "DensityGridSize": KISI_RAPAT,
              "SampleDensityGridEvery": 100,
              "WriteDensityGridEvery": max(5000, siklus // 5),
              "ComputeConventionalRDF": bool(grid),
              "NumberOfBinsConventionalRDF": 250,
              "RangeConventionalRDF": 12.0,
              "SampleConventionalRDFEvery": 100,
              "WriteConventionalRDFEvery": max(5000, siklus // 5)}
    json.dump({"SimulationType": "MonteCarlo",
               "NumberOfCycles": siklus,
               "NumberOfInitializationCycles": max(2000, siklus // 10),
               "NumberOfEquilibrationCycles": max(2000, siklus // 10),
               "PrintEvery": max(1000, siklus // 20),
               "ForceField": ".",
               "Systems": [sistem],
               "Components": komp},
              open(os.path.join(keluar, 'simulation.json'), 'w'), indent=1)
    return nama, n, len(tipe), x


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    rh_list = ([float(sys.argv[3])] if len(sys.argv) > 3 and sys.argv[3] != 'semua'
               else RH_DAFTAR)
    siklus = int(sys.argv[4]) if len(sys.argv) > 4 else 50000
    suhu = float(sys.argv[5]) if len(sys.argv) > 5 else 298.0
    x_co2 = float(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[6] != '-' else None
    grid = (sys.argv[7].lower() == 'grid') if len(sys.argv) > 7 else False
    manifes, lewat = [], 0
    for rh in rh_list:
        x = fraksi_rh(rh, x_co2)
        print('RH %3g %%  p_H2O %8.1f Pa  p_CO2 %9.2f Pa  %s' % (
            rh, x.get('H2O', 0.0) * P_TOTAL, x['CO2'] * P_TOTAL,
            'grid+RDF, swap baku' if grid else 'tanpa grid, CFCMC'))
        for f in sorted(os.listdir(src)):
            if not f.endswith('.cif'):
                continue
            # fractional RH values need their own label: int() would map
            # 2.5 % onto 2 % and 7.5 % onto 7 %.
            lab = ('%02d' % rh if float(rh).is_integer()
                   else ('%04.1f' % rh).replace('.', 'p'))
            tag = '%s__rh%s%s' % (f[:-4], lab,
                                  '' if x_co2 is None else '__co2%g' % (x_co2 * 100))
            try:
                nm, n, nt, xx = tulis(os.path.join(src, f), os.path.join(dst, tag),
                                      rh, siklus, suhu, x_co2=x_co2, grid=grid)
            except ValueError as e:
                print('  LEWATI %s : %s' % (tag, e)); lewat += 1; continue
            manifes.append(dict(tag=tag, nama=nm, RH=rh,
                                x_CO2=('%.6g' % xx['CO2']),
                                p_H2O_Pa='%.6g' % (xx.get('H2O', 0.0) * P_TOTAL),
                                p_CO2_Pa='%.6g' % (xx['CO2'] * P_TOTAL),
                                grid=('ya' if grid else ''),
                                superset='x'.join(map(str, n)), n_tipe=nt,
                                suhu=suhu, siklus=siklus))
    mp = os.path.join(dst, 'manifes_pindai.csv')
    ada = os.path.exists(mp)
    with io.open(mp, 'a' if ada else 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['tag', 'nama', 'RH', 'x_CO2', 'p_H2O_Pa',
                                           'p_CO2_Pa', 'grid', 'superset', 'n_tipe',
                                           'suhu', 'siklus'])
        if not ada:
            w.writeheader()
        for r in manifes:
            w.writerow(r)
    print('%d lari disiapkan, dilewati %d -> %s' % (len(manifes), lewat, dst))
