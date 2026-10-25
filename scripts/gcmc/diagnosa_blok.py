"""Diagnose each error bar: can it be reduced by running more cycles, or not?

Not all GCMC/Widom error bars behave in the same way. From the five production
blocks written by RASPA, three patterns can be distinguished, and only one of
them follows the 1/sqrt(cycles) law:

  DERAU   (noise) blocks scattered evenly, none dominating, no trend.
          The error decreases as 1/sqrt(cycles) and can be reduced at a
          predictable cost.

  LONJAK  (spike) one block far above the others (heavy tail). K_H/loading is
          set by a handful of rare events. Adding cycles does not reliably
          reduce the error: the tail is sampled more often and the value can
          jump. The current value is probably biased low.

  HANYUT  (drift) blocks increase or decrease monotonically. Not yet
          equilibrated. Adding cycles shifts the value rather than tightening it.

Meaning of the RASPA +/- value. Checked on 22 runs: the ratio of the reported
+/- to the standard error of the mean is 2.7760 on average, matching Student's t
for 4 degrees of freedom at 95 % (2.7764). The value reported by RASPA is
therefore a 95 % confidence interval, not 1 sigma; divide it by 2.776 to obtain sigma.

usage: python3 diagnosa_blok.py <dir> [output_name.csv]
"""
import os, re, sys, csv, io, glob, math, statistics as st

T95_4 = 2.7764                      # Student's t, 4 degrees of freedom, two-sided 95 %
AMBANG_LONJAK = 1.5                 # max/median above this -> heavy tail

PAT_KOMP = re.compile(r'^Component\s+(\d+)\s+\((\w+)\)\s*$', re.M)
PAT_BLOK = re.compile(r'Block\[\s*(\d+)\]\s+([-\d.eE+]+)')
# mixture/isotherm loading
PAT_MUAT = re.compile(
    r'Abs\. loading average\s+([-\d.eE+]+)\s+\+/-\s+([-\d.eE+]+)\s+\[mol/kg-framework\]')
# Widom
PAT_HENRY = re.compile(
    r'Henry coefficient based on Rosenbluth weight:(.*?)'
    r'Average Henry coefficient:\s+([\d.eE+-]+)\s+\+/-\s+([\d.eE+-]+)\s+\[mol/kg/Pa\]', re.S)


def spearman(v):
    n = len(v)
    urut = sorted(range(n), key=lambda i: v[i])
    pkt = [0] * n
    for r, i in enumerate(urut):
        pkt[i] = r
    d2 = sum((i - pkt[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n * n - 1))


def golongkan(b):
    """-> (pattern, max/median, relative sigma in percent, Spearman rho)"""
    if len(b) < 5 or min(b) <= 0:
        return None
    med = st.median(b)
    mm = max(b) / med if med > 0 else float('inf')
    rata = st.mean(b)
    sig = st.stdev(b) / math.sqrt(len(b))          # standard error, 1 sigma
    rel = 100 * sig / rata if rata > 0 else float('inf')
    rho = spearman(b)
    if mm > AMBANG_LONJAK:
        return 'LONJAK', mm, rel, rho
    if abs(rho) == 1.0:
        return 'HANYUT', mm, rel, rho
    return 'DERAU', mm, rel, rho


def blok_sebelum(teks, pos, n=5):
    b = [(int(m.group(1)), float(m.group(2)))
         for m in PAT_BLOK.finditer(teks, 0, pos)]
    if len(b) < n:
        return None
    lima = b[-n:]
    if [i for i, _ in lima] != list(range(n)):
        return None
    return [v for _, v in lima]


def satu_berkas(path):
    """-> [(quantity, gas, value, ci95, blocks)]"""
    teks = open(path, encoding='utf-8', errors='replace').read()
    out = []
    m = PAT_HENRY.search(teks)
    if m:
        b = [float(x) for x in re.findall(r'Block\[\s*\d+\]\s+([-\d.eE+]+)', m.group(1))]
        if len(b) >= 5:
            out.append(('K_H', '', float(m.group(2)), float(m.group(3)), b[:5]))
    judul = [(x.start(), x.group(2)) for x in PAT_KOMP.finditer(teks)]
    for x in PAT_MUAT.finditer(teks):
        gas = None
        for pos, g in judul:
            if pos < x.start():
                gas = g
            else:
                break
        b = blok_sebelum(teks, x.start())
        if b:
            out.append(('muatan', gas or '?', float(x.group(1)), float(x.group(2)), b))
    return out


if __name__ == '__main__':
    akar = os.path.expanduser(sys.argv[1])
    keluar = sys.argv[2] if len(sys.argv) > 2 else None
    baris = []
    for f in sorted(glob.glob(os.path.join(akar, '*', 'output', '*.txt'))):
        tag = os.path.basename(os.path.dirname(os.path.dirname(f)))
        try:
            hasil = satu_berkas(f)
        except Exception:
            continue
        for besaran, gas, nilai, ci, b in hasil:
            g = golongkan(b)
            if not g:
                continue
            pola, mm, rel, rho = g
            baris.append(dict(tag=tag, besaran=besaran, gas=gas,
                              nilai='%.6e' % nilai,
                              ci95_pct='%.2f' % (100 * ci / nilai) if nilai > 0 else '',
                              sigma_pct='%.2f' % rel,
                              maks_per_median='%.2f' % mm,
                              rho='%+.1f' % rho, pola=pola))
    if not baris:
        raise SystemExit('tidak ada yang terbaca dari %s' % akar)
    n = len(baris)
    print('=' * 84)
    print('DIAGNOSA BLOK — %s   (%d pengamatan)' % (akar, n))
    print('=' * 84)
    for p in ('DERAU', 'LONJAK', 'HANYUT'):
        s = [b for b in baris if b['pola'] == p]
        if not s:
            continue
        sg = sorted(float(b['sigma_pct']) for b in s)
        print('%-8s %5d (%4.1f %%)   sigma: median %5.2f %%  p90 %6.2f %%'
              % (p, len(s), 100 * len(s) / n, sg[len(sg) // 2], sg[int(.9 * (len(sg) - 1))]))
    print()
    print('  DERAU  bisa ditekan: sigma turun 1/sqrt(siklus)')
    print('  LONJAK tidak bisa diandalkan; nilainya mungkin bias rendah')
    print('  HANYUT belum setimbang; siklus tambahan MENGGESER nilai')
    if keluar:
        with io.open(keluar, 'w', encoding='utf-8', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(baris[0].keys()))
            w.writeheader()
            for b in baris:
                w.writerow(b)
        print('\n-> %s' % keluar)
