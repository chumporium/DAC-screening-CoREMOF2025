"""Framework heat capacity from the CoRE MOF 2025 metadata, for all 7,820 structures.

Two sources:
  - CSD-modified : Heat_capacity@{300,350,400}K (J/g/K) + std, plus cp0 and k_cp
  - SI subset    : heat_capacity {300K, 350K, 400K}, without std and without a fit

The values are predictions of an ML model (RAC features), not measurements. They
are useful for ranking, not as absolute numbers.

The output stores cp0 and k_cp so that Cp(T) = cp0 + k_cp x T can be integrated
between the adsorption and regeneration temperatures. For SI entries without a fit,
the line is fitted here to the three points (least squares).
"""
import os, csv, io, json, math

D = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(D, 'data')
CSD = os.environ.get('COREMOF_CSD_CSV', 'CR_data_CSD_modified_20250227.csv')      # CoRE MOF 2025 metadata table
SI_JSON = os.path.join(DATA, 'meta_si.json')
OUT = os.path.join(DATA, 'cp_semua.csv')

KOL = ['nama', 'Cp_300K', 'Cp_std_300K', 'Cp_350K', 'Cp_400K', 'cp0', 'k_cp', 'sumber_Cp']
T = (300.0, 350.0, 400.0)


def fl(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def tanpa_pacman(k):
    return k[:-7] if k.endswith('_pacman') else k


def garis(ys):
    """Least-squares line through the three points (300/350/400 K) -> (cp0, k_cp)."""
    n = len(ys)
    mx, my = sum(T) / n, sum(ys) / n
    sxx = sum((t - mx) ** 2 for t in T)
    if sxx == 0:
        return my, 0.0
    k = sum((t - mx) * (y - my) for t, y in zip(T, ys)) / sxx
    return my - k * mx, k


def kumpulkan():
    d = {}
    if os.path.exists(CSD):
        for r in csv.DictReader(io.open(CSD, encoding='utf-8-sig')):
            c = [fl(r.get('Heat_capacity@%dK (J/g/K)' % t)) for t in (300, 350, 400)]
            if c[0] is None:
                continue
            cp0, kcp = fl(r.get('cp0 (J/g/K)')), fl(r.get('k_cp (J/g/K/K)'))
            if cp0 is None or kcp is None:
                cp0, kcp = garis([x for x in c if x is not None] or [c[0]])
            d[tanpa_pacman(r['refcode'])] = dict(
                Cp_300K=c[0], Cp_std_300K=fl(r.get('std @ 300 K (J/g/K)')),
                Cp_350K=c[1], Cp_400K=c[2], cp0=cp0, k_cp=kcp,
                sumber_Cp='CSD-modified')
    if os.path.exists(SI_JSON):
        for k, v in json.load(open(SI_JSON, encoding='utf-8')).items():
            k = tanpa_pacman(k)
            if k in d:
                continue
            h = v.get('heat_capacity')
            if not isinstance(h, dict) or '300K' not in h:
                continue
            c = [fl(h.get('%dK' % t)) for t in (300, 350, 400)]
            if c[0] is None:
                continue
            cp0, kcp = garis([x for x in c if x is not None])
            d[k] = dict(Cp_300K=c[0], Cp_std_300K=None, Cp_350K=c[1], Cp_400K=c[2],
                        cp0=cp0, k_cp=kcp, sumber_Cp='SI CoRE MOF')
    return d


if __name__ == '__main__':
    CP = kumpulkan()
    nama = [r['nama'] for r in csv.DictReader(
        open(os.path.join(DATA, 'hasil_final.csv'), encoding='utf-8'))]
    with io.open(OUT, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=KOL)
        w.writeheader()
        n = 0
        for x in nama:
            c = CP.get(x)
            if not c:
                continue
            w.writerow(dict(nama=x, **c))
            n += 1
    v = sorted(c['Cp_300K'] for c in CP.values())
    print('-> %s' % OUT)
    print('   dalam metadata      : %d struktur' % len(CP))
    print('   cocok ke basis kita : %d dari %d  (%.1f %%)' % (n, len(nama), 100 * n / len(nama)))
    print('   Cp@300 K            : %.3f - %.3f, median %.3f J/(g K)'
          % (v[0], v[-1], v[len(v) // 2]))
