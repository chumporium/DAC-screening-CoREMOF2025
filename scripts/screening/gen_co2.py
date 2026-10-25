"""RASPA3 Widom input generator for the K_H of CO2, with a selectable CO2 model.

usage: python3 gen_co2.py <dir_cif> <dir_out> <model> [cycles]
model: trappe | garcia
"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_widom import UFF, MASSA, baca_cif, lebar_tegak

CO2 = {
    # sigma [A], epsilon/kB [K], charge [e], C=O bond length [A]
    'trappe': dict(sigC=2.80,  epsC=27.0,    qC=0.70,   sigO=3.05,  epsO=79.0,   qO=-0.35,   l=1.16,
                   sumber='TraPPE, Potoff & Siepmann, AIChE J 47, 1676 (2001)'),
    'garcia': dict(sigC=2.745, epsC=29.933,  qC=0.6512, sigO=3.017, epsO=85.671, qO=-0.3256, l=1.149,
                   sumber='Garcia-Sanchez et al., J Phys Chem C 113, 8814 (2009)'),
}


def tulis(cif, keluar, model, cutoff=9.0, siklus=2000, suhu=298.0):
    m = CO2[model]
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

    if not atom:
        raise ValueError('CIF tanpa atom terbaca: %s' % cif)
    hilang = sorted({a['el'] for a in atom} - set(UFF))
    if hilang:
        raise ValueError('unsur tanpa parameter UFF: %s' % ','.join(hilang))
    ps, si = [], []
    for i, a in enumerate(atom):
        el = a['el']
        ps.append({"name": "A%d" % i, "framework": True, "print_to_output": False, "element": el,
                   "print_as": el, "mass": MASSA.get(el, 1.0), "charge": a['q'], "source": "PACMAN-DDEC6"})
        si.append({"name": "A%d" % i, "type": "lennard-jones", "parameters": list(UFF[el]), "source": "UFF"})
    ps += [
        {"name": "C_co2", "framework": False, "print_to_output": True, "element": "C", "print_as": "C",
         "mass": 12.0, "charge": m['qC'], "source": m['sumber']},
        {"name": "O_co2", "framework": False, "print_to_output": True, "element": "O", "print_as": "O",
         "mass": 15.9994, "charge": m['qO'], "source": m['sumber']}]
    si += [
        {"name": "C_co2", "type": "lennard-jones", "parameters": [m['epsC'], m['sigC']], "source": m['sumber']},
        {"name": "O_co2", "type": "lennard-jones", "parameters": [m['epsO'], m['sigO']], "source": m['sumber']}]

    json.dump({"MixingRule": "Lorentz-Berthelot", "TruncationMethod": "shifted",
               "TailCorrections": False, "PseudoAtoms": ps, "SelfInteractions": si},
              open(os.path.join(keluar, 'force_field.json'), 'w'), indent=1)
    json.dump({"CriticalTemperature": 304.1282, "CriticalPressure": 7377300.0,
               "AcentricFactor": 0.22394, "Type": "rigid",
               "pseudoAtoms": [["O_co2", [0.0, 0.0, m['l']]], ["C_co2", [0.0, 0.0, 0.0]],
                               ["O_co2", [0.0, 0.0, -m['l']]]],
               "Bonds": [[0, 1], [1, 2]]},
              open(os.path.join(keluar, 'CO2.json'), 'w'), indent=1)

    n = [max(1, int(math.ceil(2 * cutoff / w))) for w in lebar_tegak(sel)]
    json.dump({"SimulationType": "MonteCarlo", "NumberOfCycles": siklus,
               "NumberOfInitializationCycles": 0, "PrintEvery": max(1, siklus // 4),
               "ForceField": ".",
               "Systems": [{"Type": "Framework", "Name": nama, "NumberOfUnitCells": n,
                            "ExternalTemperature": suhu, "ExternalPressure": 1.0e5,
                            "ChargeMethod": "Ewald", "CutOffVDW": cutoff, "CutOffCoulomb": cutoff}],
               "Components": [{"Name": "CO2", "WidomProbability": 1.0, "CreateNumberOfMolecules": 0}]},
              open(os.path.join(keluar, 'simulation.json'), 'w'), indent=1)
    return nama, n


if __name__ == '__main__':
    src, dst, model = sys.argv[1], sys.argv[2], sys.argv[3]
    siklus = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    dilewati = []
    for f in sorted(os.listdir(src)):
        if not f.endswith('.cif'):
            continue
        try:
            nm, n = tulis(os.path.join(src, f), os.path.join(dst, f[:-4]), model, siklus=siklus)
        except ValueError as e:
            print('LEWATI %-34s %s' % (f[:-4], e)); dilewati.append((f[:-4], str(e))); continue
        print('%-8s %-34s superset=%s' % (model, nm, n))
    if dilewati:
        print('TOTAL DILEWATI: %d' % len(dilewati))
