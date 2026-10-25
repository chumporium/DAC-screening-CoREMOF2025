import sys, os, re, json, math

UFF = {  # epsilon/kB [K], sigma [A]  (Rappe 1992)
 'H':(22.142,2.571),'C':(52.838,3.431),'N':(34.722,3.261),'O':(30.193,3.118),
 'F':(25.160,2.997),'Si':(202.294,3.826),'S':(137.882,3.595),'Cl':(114.231,3.516),
 'Sc':(9.561,2.936),'Ti':(8.554,2.829),'V':(8.051,2.801),'Cr':(7.548,2.693),
 'Mn':(6.541,2.638),'Fe':(6.541,2.594),'Co':(7.045,2.559),'Ni':(7.548,2.525),
 'Cu':(2.516,3.114),'Zn':(62.399,2.462),'Zr':(34.722,2.783),'Nb':(29.710,2.820),
 'Cd':(114.720,2.537),'Ag':(18.116,2.805),'Br':(126.3,3.732),'Al':(155.998,4.008),
 'Mg':(55.857,2.691),'Ca':(119.766,3.028),'Na':(15.096,2.658),'K':(17.612,3.396),
 'In':(301.428,3.976),'Mo':(28.203,2.719),'Ru':(28.203,2.640),'Rh':(26.688,2.609),
 'Pd':(24.173,2.583),'Pt':(40.257,2.454),'Au':(19.632,2.934),'Pb':(333.634,3.828),
 'Tb':(3.523,3.243),'Eu':(4.028,3.112),'La':(8.554,3.138),'Ce':(6.541,3.168),
 'Nd':(5.032,3.185),'Sm':(4.028,3.136),'Gd':(4.531,3.000),'Dy':(3.523,3.054),
 'Ho':(3.523,3.037),'Er':(3.523,3.021),'Yb':(114.72,2.989),'Y':(36.229,2.980),
 'Sr':(118.258,3.244),'Ba':(183.185,3.299),'Li':(12.581,2.184),'Be':(42.774,2.446),
 'B':(90.580,3.638),'P':(153.482,3.695),'As':(155.495,3.769),'Se':(146.437,3.746),
 'Sb':(225.946,3.938),'Te':(200.281,3.982),'I':(170.591,4.009),'Ge':(190.72,3.813),
 'Sn':(285.34,3.913),'Hf':(36.229,2.798),'Ta':(40.760,2.824),'W':(33.716,2.734),
 'Re':(33.213,2.632),'Os':(18.619,2.780),'Ir':(36.732,2.530),'Hg':(193.740,2.410),
 'Bi':(260.669,3.893),'U':(11.071,3.225),'Th':(13.083,3.383),
 'Ga':(208.836,3.905),'Tl':(342.192,3.873),'Rb':(20.128,3.665),'Cs':(22.644,4.024),
 'Tc':(24.173,2.671),'Pr':(5.032,3.213),'Pm':(4.531,3.160),'Tm':(3.019,3.006),
 'Lu':(20.630,3.243),'Np':(9.561,3.050),'Pu':(8.051,3.050),
}
MASSA = {'H':1.008,'C':12.011,'N':14.007,'O':15.999,'F':18.998,'Si':28.086,'S':32.06,
 'Cl':35.45,'Sc':44.956,'Ti':47.867,'V':50.942,'Cr':51.996,'Mn':54.938,'Fe':55.845,
 'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,'Zr':91.224,'Nb':92.906,'Cd':112.411,
 'Ag':107.868,'Br':79.904,'Al':26.982,'Mg':24.305,'Ca':40.078,'Na':22.990,'K':39.098,
 'Y':88.906,'Mo':95.95,'La':138.905,'Ce':140.116,'Nd':144.242,'Sm':150.36,
 'Eu':151.964,'Gd':157.25,'Tb':158.925,'Dy':162.5,'Ho':164.930,'Er':167.259,
 'Yb':173.045,'Ba':137.327,'Sr':87.62,'Li':6.94,'Be':9.012,'B':10.81,'P':30.974,
 'As':74.922,'Se':78.971,'Sb':121.760,'Te':127.60,'I':126.904,'Ge':72.63,'Sn':118.71,
 'Hf':178.49,'Ta':180.948,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,
 'Pt':195.084,'Au':196.967,'Hg':200.592,'Pb':207.2,'Bi':208.980,'U':238.029,
 'Th':232.038,'In':114.818,'Ru':101.07,'Rh':102.906,'Pd':106.42,
 'Ga':69.723,'Tl':204.38,'Rb':85.468,'Cs':132.905,'Tc':98.0,'Pr':140.908,
 'Pm':145.0,'Tm':168.934,'Lu':174.967,'Np':237.0,'Pu':244.0}


def baca_cif(p):
    txt = open(p, encoding='utf-8', errors='ignore').read()

    def num(tag):
        m = re.search(r'^\s*' + tag + r'\s+([\d.eE+-]+)', txt, re.M)
        return float(m.group(1))

    sel = [num('_cell_length_a'), num('_cell_length_b'), num('_cell_length_c'),
           num('_cell_angle_alpha'), num('_cell_angle_beta'), num('_cell_angle_gamma')]
    atom = []
    for l in txt.splitlines():
        s = l.split()
        if len(s) == 8 and re.match(r'^[A-Z][a-z]?$', s[0]) and re.match(r'^[A-Z][a-z]?\d', s[1]):
            atom.append({'el': s[0], 'x': s[3], 'y': s[4], 'z': s[5], 'occ': s[6], 'q': float(s[7])})
    return sel, atom


def lebar_tegak(sel):
    a, b, c, al, be, ga = sel
    al, be, ga = [math.radians(v) for v in (al, be, ga)]
    V = a * b * c * math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2
                              + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    return [V / (b * c * math.sin(al)), V / (a * c * math.sin(be)), V / (a * b * math.sin(ga))]


def buat(cif, keluar, cutoff=9.0, siklus=2000, suhu=298.0):
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
        if el not in UFF:
            raise SystemExit('UFF tidak punya %s (%s)' % (el, nama))
        ps.append({"name": "A%d" % i, "framework": True, "print_to_output": False, "element": el,
                   "print_as": el, "mass": MASSA.get(el, 1.0), "charge": a['q'], "source": "PACMAN-DDEC6"})
        si.append({"name": "A%d" % i, "type": "lennard-jones", "parameters": list(UFF[el]), "source": "UFF"})
    ps += [
        {"name": "Ow", "framework": False, "print_to_output": True, "element": "O", "print_as": "O",
         "mass": 15.9994, "charge": 0.0, "source": "TIP4P"},
        {"name": "Hw", "framework": False, "print_to_output": True, "element": "H", "print_as": "H",
         "mass": 1.00794, "charge": 0.52, "source": "TIP4P"},
        {"name": "Mw", "framework": False, "print_to_output": False, "element": "O", "print_as": "-",
         "mass": 0.0, "charge": -1.04, "source": "TIP4P"}]
    si += [
        {"name": "Ow", "type": "lennard-jones", "parameters": [78.0, 3.15365], "source": "TIP4P"},
        {"name": "Hw", "type": "none", "parameters": [0.0, 1.0], "source": "TIP4P"},
        {"name": "Mw", "type": "none", "parameters": [0.0, 1.0], "source": "TIP4P"}]
    json.dump({"MixingRule": "Lorentz-Berthelot", "TruncationMethod": "shifted",
               "TailCorrections": False, "PseudoAtoms": ps, "SelfInteractions": si},
              open(os.path.join(keluar, 'force_field.json'), 'w'), indent=1)

    json.dump({"CriticalTemperature": 647.14, "CriticalPressure": 22064000.0,
               "AcentricFactor": 0.3449, "Type": "rigid",
               "pseudoAtoms": [["Ow", [0.0, 0.0, 0.0]], ["Hw", [0.7569, 0.5859, 0.0]],
                               ["Hw", [-0.7569, 0.5859, 0.0]], ["Mw", [0.0, 0.15, 0.0]]],
               "Bonds": [[0, 1], [0, 2], [0, 3]]},
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
    src, dst = sys.argv[1], sys.argv[2]
    siklus = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    for f in sorted(os.listdir(src)):
        if not f.endswith('.cif'):
            continue
        nm, n, na = buat(os.path.join(src, f), os.path.join(dst, f[:-4]), siklus=siklus)
        print('%-34s superset=%s atom_sel=%d total=%d' % (nm, n, na, na * n[0] * n[1] * n[2]))
