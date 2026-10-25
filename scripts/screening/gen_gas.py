"""RASPA3 Widom input generator for N2 / O2 / H2O / CO2 / Ar in CoRE MOF frameworks.

Framework: UFF LJ + per-atom PACMAN-DDEC6 charges from the CIF.
usage: python3 gen_gas.py <dir_cif> <dir_out> <gas> [cycles] [T_K]
gas: co2 | n2 | o2 | h2o | ar
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_widom import UFF, MASSA, baca_cif, lebar_tegak

# sites: (name, sigma A, eps/kB K, charge e, z A, mass)  -- LJ type 'none' if sigma is None
GAS = {
    'co2': dict(Tc=304.1282, Pc=7377300.0, w=0.22394, sumber='TraPPE (Potoff & Siepmann 2001)',
                situs=[('O_co2', 3.05, 79.0, -0.35, +1.16, 15.9994),
                       ('C_co2', 2.80, 27.0, +0.70, 0.0, 12.0),
                       ('O_co2', 3.05, 79.0, -0.35, -1.16, 15.9994)],
                ikatan=[[0, 1], [1, 2]]),
    'n2': dict(Tc=126.192, Pc=3395800.0, w=0.0372, sumber='TraPPE (Potoff & Siepmann 2001)',
               situs=[('N_n2', 3.31, 36.0, -0.482, +0.55, 14.00674),
                      ('N_com', None, None, +0.964, 0.0, 0.0),
                      ('N_n2', 3.31, 36.0, -0.482, -0.55, 14.00674)],
               ikatan=[[0, 1], [1, 2]]),
    'o2': dict(Tc=154.581, Pc=5043000.0, w=0.0222, sumber='TraPPE (Potoff & Siepmann 2001)',
               situs=[('O_o2', 3.013, 49.1, -0.123, +0.605, 15.9994),
                      ('O_com', None, None, +0.246, 0.0, 0.0),
                      ('O_o2', 3.013, 49.1, -0.123, -0.605, 15.9994)],
               ikatan=[[0, 1], [1, 2]]),
    'ar': dict(Tc=150.687, Pc=4863000.0, w=-0.002, sumber='LJ Ar (sigma 3,405 A; eps/kB 119,8 K)',
               situs=[('Ar', 3.405, 119.8, 0.0, 0.0, 39.948)],
               ikatan=[]),
    'h2o': dict(Tc=647.14, Pc=22064000.0, w=0.3449, sumber='TIP4P',
                situs=[('Ow', 3.15365, 78.0, 0.0, None, 15.9994),
                       ('Hw', None, None, +0.52, None, 1.00794),
                       ('Hw', None, None, +0.52, None, 1.00794),
                       ('Mw', None, None, -1.04, None, 0.0)],
                ikatan=[[0, 1], [0, 2], [0, 3]],
                koord=[[0.0, 0.0, 0.0], [0.7569, 0.5859, 0.0],
                       [-0.7569, 0.5859, 0.0], [0.0, 0.15, 0.0]]),
}


def tulis(cif, keluar, gas, cutoff=9.0, siklus=2000, suhu=298.0):
    m = GAS[gas]
    nama = os.path.basename(cif)[:-4]
    sel, atom = baca_cif(cif)
    if not atom:
        raise ValueError('CIF tanpa atom: %s' % cif)
    hilang = sorted({a['el'] for a in atom} - set(UFF))
    if hilang:
        raise ValueError('tanpa UFF: ' + ','.join(hilang))
    os.makedirs(keluar, exist_ok=True)
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
        for i, a in enumerate(atom):
            f.write('%-3s A%-6d 1.0  %s %s %s %s\n' % (a['el'], i, a['x'], a['y'], a['z'], a['occ']))

    ps, si = [], []
    for i, a in enumerate(atom):
        el = a['el']
        ps.append({"name": "A%d" % i, "framework": True, "print_to_output": False, "element": el,
                   "print_as": el, "mass": MASSA.get(el, 1.0), "charge": a['q'],
                   "source": "PACMAN-DDEC6"})
        si.append({"name": "A%d" % i, "type": "lennard-jones", "parameters": list(UFF[el]),
                   "source": "UFF"})
    sudah = set()
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

    if 'koord' in m:
        situs = [[s[0], k] for s, k in zip(m['situs'], m['koord'])]
    else:
        situs = [[s[0], [0.0, 0.0, s[4]]] for s in m['situs']]
    json.dump({"CriticalTemperature": m['Tc'], "CriticalPressure": m['Pc'],
               "AcentricFactor": m['w'], "Type": "rigid",
               "pseudoAtoms": situs, "Bonds": m['ikatan']},
              open(os.path.join(keluar, gas.upper() + '.json'), 'w'), indent=1)

    n = [max(1, int(math.ceil(2 * cutoff / w))) for w in lebar_tegak(sel)]
    json.dump({"SimulationType": "MonteCarlo", "NumberOfCycles": siklus,
               "NumberOfInitializationCycles": 0, "PrintEvery": max(1, siklus // 4),
               "ForceField": ".",
               "Systems": [{"Type": "Framework", "Name": nama, "NumberOfUnitCells": n,
                            "ExternalTemperature": suhu, "ExternalPressure": 1.0e5,
                            "ChargeMethod": "Ewald", "CutOffVDW": cutoff,
                            "CutOffCoulomb": cutoff}],
               "Components": [{"Name": gas.upper(), "WidomProbability": 1.0,
                               "CreateNumberOfMolecules": 0}]},
              open(os.path.join(keluar, 'simulation.json'), 'w'), indent=1)
    return nama, n


if __name__ == '__main__':
    src, dst, gas = sys.argv[1], sys.argv[2], sys.argv[3].lower()
    siklus = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    suhu = float(sys.argv[5]) if len(sys.argv) > 5 else 298.0
    lewat = 0
    for f in sorted(os.listdir(src)):
        if not f.endswith('.cif'):
            continue
        try:
            tulis(os.path.join(src, f), os.path.join(dst, f[:-4]), gas, siklus=siklus, suhu=suhu)
        except ValueError as e:
            print('LEWATI %s : %s' % (f[:-4], e))
            lewat += 1
    print('%s selesai, dilewati %d' % (gas, lewat))
