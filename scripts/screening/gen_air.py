"""RASPA3 Widom input generator for the K_H of water, for several water models.

usage: python3 gen_air.py <dir_cif> <dir_out> <model> [cycles]
model: tip4p | tip4p2005 | spce | tip3p
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_widom import UFF, MASSA, baca_cif, lebar_tegak

# geometry: each H at half the H-O-H angle from the y axis
def h_pos(r_oh, sudut):
    a = math.radians(sudut / 2.0)
    return [round(r_oh * math.sin(a), 5), round(r_oh * math.cos(a), 5), 0.0]

AIR = {
    'tip4p':     dict(eps=78.0,  sig=3.15365, qO=0.0,     qH=0.52,   qM=-1.04,   dM=0.15,   r=0.9572, ang=104.52),
    'tip4p2005': dict(eps=93.2,  sig=3.1589,  qO=0.0,     qH=0.5564, qM=-1.1128, dM=0.1546, r=0.9572, ang=104.52),
    'spce':      dict(eps=78.20, sig=3.166,   qO=-0.8476, qH=0.4238, qM=None,    dM=None,   r=1.0,    ang=109.47),
    'tip3p':     dict(eps=76.54, sig=3.1506,  qO=-0.834,  qH=0.417,  qM=None,    dM=None,   r=0.9572, ang=104.52),
}


def tulis(cif, keluar, model, cutoff=9.0, siklus=2000, suhu=298.0):
    m = AIR[model]
    nama = os.path.basename(cif)[:-4]
    sel, atom = baca_cif(cif)
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
                   "print_as": el, "mass": MASSA.get(el, 1.0), "charge": a['q'], "source": "PACMAN-DDEC6"})
        si.append({"name": "A%d" % i, "type": "lennard-jones", "parameters": list(UFF[el]), "source": "UFF"})

    ps.append({"name": "Ow", "framework": False, "print_to_output": True, "element": "O",
               "print_as": "O", "mass": 15.9994, "charge": m['qO'], "source": model})
    ps.append({"name": "Hw", "framework": False, "print_to_output": True, "element": "H",
               "print_as": "H", "mass": 1.00794, "charge": m['qH'], "source": model})
    si.append({"name": "Ow", "type": "lennard-jones", "parameters": [m['eps'], m['sig']], "source": model})
    si.append({"name": "Hw", "type": "none", "parameters": [0.0, 1.0], "source": model})

    situs = [["Ow", [0.0, 0.0, 0.0]],
             ["Hw", h_pos(m['r'], m['ang'])],
             ["Hw", [-h_pos(m['r'], m['ang'])[0], h_pos(m['r'], m['ang'])[1], 0.0]]]
    ikat = [[0, 1], [0, 2]]
    if m['qM'] is not None:
        ps.append({"name": "Mw", "framework": False, "print_to_output": False, "element": "O",
                   "print_as": "-", "mass": 0.0, "charge": m['qM'], "source": model})
        si.append({"name": "Mw", "type": "none", "parameters": [0.0, 1.0], "source": model})
        situs.append(["Mw", [0.0, m['dM'], 0.0]])
        ikat.append([0, 3])

    json.dump({"MixingRule": "Lorentz-Berthelot", "TruncationMethod": "shifted",
               "TailCorrections": False, "PseudoAtoms": ps, "SelfInteractions": si},
              open(os.path.join(keluar, 'force_field.json'), 'w'), indent=1)
    json.dump({"CriticalTemperature": 647.14, "CriticalPressure": 22064000.0,
               "AcentricFactor": 0.3449, "Type": "rigid",
               "pseudoAtoms": situs, "Bonds": ikat},
              open(os.path.join(keluar, 'water.json'), 'w'), indent=1)

    n = [max(1, int(math.ceil(2 * cutoff / w))) for w in lebar_tegak(sel)]
    json.dump({"SimulationType": "MonteCarlo", "NumberOfCycles": siklus,
               "NumberOfInitializationCycles": 0, "PrintEvery": max(1, siklus // 4),
               "ForceField": ".",
               "Systems": [{"Type": "Framework", "Name": nama, "NumberOfUnitCells": n,
                            "ExternalTemperature": suhu, "ExternalPressure": 1.0e5,
                            "ChargeMethod": "Ewald", "CutOffVDW": cutoff, "CutOffCoulomb": cutoff}],
               "Components": [{"Name": "water", "WidomProbability": 1.0, "CreateNumberOfMolecules": 0}]},
              open(os.path.join(keluar, 'simulation.json'), 'w'), indent=1)
    return nama, n, len(atom)


if __name__ == '__main__':
    src, dst, model = sys.argv[1], sys.argv[2], sys.argv[3]
    siklus = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    for f in sorted(os.listdir(src)):
        if not f.endswith('.cif'):
            continue
        nm, n, na = tulis(os.path.join(src, f), os.path.join(dst, f[:-4]), model, siklus=siklus)
        print('%-10s %-34s superset=%s' % (model, nm, n))
