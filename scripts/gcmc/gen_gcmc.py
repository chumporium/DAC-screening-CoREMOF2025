"""RASPA3 GCMC input generator: realistic air mixture at 1 atm.

Framework: UFF LJ + PACMAN-DDEC6 charges; pseudo-atoms are grouped by
(element, charge rounded to 3 decimals).

usage: python3 gen_gcmc.py <dir_cif> <dir_out> <condition> [cycles] [T_K] [mode]
  condition : kering (dry) | basah (humid)
  mode      : baku (standard) | cfcmc   (cfcmc uses CFCMC_CBMC_SwapProbability for CO2)

Why CFCMC: at 1 atm CO2 is only 0.04 % of the mixture, so its partial fugacity
is 41 Pa. Direct insertions are rarely accepted and the CO2 uptake statistics
evolve very slowly, while the output still reports numbers that look well
behaved. CFCMC inserts molecules gradually through a lambda parameter, which
greatly increases the acceptance probability. The 'baku' mode is provided for
comparison: if both modes agree, CFCMC is not required.
"""
import sys, os, json, math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_widom import UFF, MASSA, baca_cif, lebar_tegak

# sites: (name, sigma A, eps/kB K, charge e, z A, mass); sigma None -> no LJ
GAS = {
    'CO2': dict(Tc=304.1282, Pc=7377300.0, w=0.22394, sumber='TraPPE',
                situs=[('O_co2', 3.05, 79.0, -0.35, +1.16, 15.9994),
                       ('C_co2', 2.80, 27.0, +0.70, 0.0, 12.0),
                       ('O_co2', 3.05, 79.0, -0.35, -1.16, 15.9994)],
                ikatan=[[0, 1], [1, 2]]),
    'N2': dict(Tc=126.192, Pc=3395800.0, w=0.0372, sumber='TraPPE',
               situs=[('N_n2', 3.31, 36.0, -0.482, +0.55, 14.00674),
                      ('N_com', None, None, +0.964, 0.0, 0.0),
                      ('N_n2', 3.31, 36.0, -0.482, -0.55, 14.00674)],
               ikatan=[[0, 1], [1, 2]]),
    'O2': dict(Tc=154.581, Pc=5043000.0, w=0.0222, sumber='TraPPE',
               situs=[('O_o2', 3.013, 49.1, -0.123, +0.605, 15.9994),
                      ('O_com', None, None, +0.246, 0.0, 0.0),
                      ('O_o2', 3.013, 49.1, -0.123, -0.605, 15.9994)],
               ikatan=[[0, 1], [1, 2]]),
    'AR': dict(Tc=150.687, Pc=4863000.0, w=-0.002, sumber='LJ Ar',
               situs=[('Ar', 3.405, 119.8, 0.0, 0.0, 39.948)], ikatan=[]),
    'H2O': dict(Tc=647.14, Pc=22064000.0, w=0.3449, sumber='TIP4P',
                situs=[('Ow', 3.15365, 78.0, 0.0, None, 15.9994),
                       ('Hw', None, None, +0.52, None, 1.00794),
                       ('Hw', None, None, +0.52, None, 1.00794),
                       ('Mw', None, None, -1.04, None, 0.0)],
                ikatan=[[0, 1], [0, 2], [0, 3]],
                koord=[[0.0, 0.0, 0.0], [0.7569, 0.5859, 0.0],
                       [-0.7569, 0.5859, 0.0], [0.0, 0.15, 0.0]]),
}

P_TOTAL = 101325.0
X_KERING = {'N2': 0.78084, 'O2': 0.20946, 'AR': 0.00934, 'CO2': 0.000415}
P_H2O_RH50 = 1585.0


def fraksi(kondisi):
    if kondisi == 'kering':
        return dict(X_KERING)
    x_air = P_H2O_RH50 / P_TOTAL
    sisa = 1.0 - x_air
    d = {k: v * sisa for k, v in X_KERING.items()}
    d['H2O'] = x_air
    return d


def tulis(cif, keluar, kondisi, siklus=50000, suhu=298.0, mode='cfcmc', cutoff=9.0):
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
    x = fraksi(kondisi)
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
        if 'koord' in m:
            situs = [[s[0], k] for s, k in zip(m['situs'], m['koord'])]
        else:
            situs = [[s[0], [0.0, 0.0, s[4]]] for s in m['situs']]
        json.dump({"CriticalTemperature": m['Tc'], "CriticalPressure": m['Pc'],
                   "AcentricFactor": m['w'], "Type": "rigid",
                   "pseudoAtoms": situs, "Bonds": m['ikatan']},
                  open(os.path.join(keluar, gas + '.json'), 'w'), indent=1)

    komp = []
    for gas, xf in sorted(x.items(), key=lambda z: -z[1]):
        c = {"Name": gas, "MolFraction": round(xf, 8), "FugacityCoefficient": 1.0,
             "TranslationProbability": 0.5, "ReinsertionProbability": 0.5,
             "WidomProbability": 0.2, "CreateNumberOfMolecules": 0}
        if gas != 'AR':
            c["RotationProbability"] = 0.5
        # CO2 uses CFCMC because its partial fugacity is only 41 Pa
        if mode == 'cfcmc' and gas == 'CO2':
            c["CFCMC_CBMC_SwapProbability"] = 1.0
            c["ThermodynamicIntegration"] = True
        else:
            c["SwapProbability"] = 1.0
        komp.append(c)

    n = [max(1, int(math.ceil(2 * cutoff / w))) for w in lebar_tegak(sel)]
    json.dump({"SimulationType": "MonteCarlo",
               "NumberOfCycles": siklus,
               "NumberOfInitializationCycles": max(2000, siklus // 10),
               "NumberOfEquilibrationCycles": max(2000, siklus // 10),
               "PrintEvery": max(1000, siklus // 20),
               "ForceField": ".",
               "Systems": [{"Type": "Framework", "Name": nama, "NumberOfUnitCells": n,
                            "ExternalTemperature": suhu, "ExternalPressure": P_TOTAL,
                            "ChargeMethod": "Ewald", "CutOffVDW": cutoff,
                            "CutOffCoulomb": cutoff}],
               "Components": komp},
              open(os.path.join(keluar, 'simulation.json'), 'w'), indent=1)
    return nama, n, len(tipe)


if __name__ == '__main__':
    src, dst, kondisi = sys.argv[1], sys.argv[2], sys.argv[3].lower()
    siklus = int(sys.argv[4]) if len(sys.argv) > 4 else 50000
    suhu = float(sys.argv[5]) if len(sys.argv) > 5 else 298.0
    mode = sys.argv[6].lower() if len(sys.argv) > 6 else 'cfcmc'
    x = fraksi(kondisi)
    print('kondisi %s pada %.0f Pa: %s' % (
        kondisi, P_TOTAL, ', '.join('%s %.6f' % (k, v) for k, v in
                                    sorted(x.items(), key=lambda z: -z[1]))))
    lewat = 0
    for f in sorted(os.listdir(src)):
        if not f.endswith('.cif'):
            continue
        try:
            nm, n, nt = tulis(os.path.join(src, f), os.path.join(dst, f[:-4]),
                              kondisi, siklus, suhu, mode)
            print('  %-24s superset %s  tipe %d' % (nm, 'x'.join(map(str, n)), nt))
        except ValueError as e:
            print('  LEWATI %s : %s' % (f[:-4], e))
            lewat += 1
    print('%s/%s selesai, dilewati %d' % (kondisi, mode, lewat))
