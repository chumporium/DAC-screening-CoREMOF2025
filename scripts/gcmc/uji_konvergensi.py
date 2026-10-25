"""Is a GCMC error bar noise (converged) or drift (insufficient production)?

The 1/sqrt(loading) scaling argument is not sufficient: the error also decreases
as 1/sqrt(cycles), so varying the loading at a fixed run length cannot show
that the run length is adequate.

The proper test uses the 5 production blocks written by RASPA:

  HANYUT  (drift) early blocks differ systematically from late blocks -> not yet
          equilibrated, production too short. The mean is biased.
  DERAU   (noise) blocks scattered without a trend -> equilibrated; adding cycles
          only shrinks the error bar and does not shift the value.

Test statistic: Spearman rank correlation between block index and block value.
With 5 blocks, |rho| = 1.0 (perfect ordering) occurs with probability
2/120 = 1.7 % for random data. Hence |rho| = 1 is classified as drift and
everything else as noise.

usage: python3 uji_konvergensi.py <run_dir> [max_runs]
"""
import os, re, sys, csv, io, statistics as st

PAT_KOMP = re.compile(r'^Component\s+(\d+)\s+\((\w+)\)\s*$', re.M)
PAT_BLOK = re.compile(r'^\s*Block\[\s*(\d+)\]\s+([-\d.eE+]+)\s*$', re.M)
PAT_ABS = re.compile(
    r'Abs\. loading average\s+([-\d.eE+]+)\s+\+/-\s+([-\d.eE+]+)\s+\[mol/kg-framework\]')


def blok_akhir(path):
    """-> {gas: ([5 block values], q, error)} from the production mol/kg-framework blocks."""
    teks = open(path, encoding='utf-8', errors='replace').read()
    judul = [(m.start(), m.group(2)) for m in PAT_KOMP.finditer(teks)]
    out = {}
    for m in PAT_ABS.finditer(teks):
        if 'mol/kg' not in m.group(0):
            continue
        gas = None
        for pos, g in judul:
            if pos < m.start():
                gas = g
            else:
                break
        if gas is None:
            continue
        # the mol/kg blocks are the 5 lines directly above the mol/kg average line;
        # take the last 5 Block[] entries before this position
        b = [(int(x.group(1)), float(x.group(2)))
             for x in PAT_BLOK.finditer(teks, 0, m.start())]
        if len(b) < 5:
            continue
        lima = b[-5:]
        if [i for i, _ in lima] != [0, 1, 2, 3, 4]:
            continue
        out[gas] = ([v for _, v in lima], float(m.group(1)), float(m.group(2)))
    return out


def spearman5(v):
    """Spearman rho between block index (0..4) and value, in [-1, 1]."""
    n = len(v)
    urut = sorted(range(n), key=lambda i: v[i])
    pkt = [0] * n
    for r, i in enumerate(urut):
        pkt[i] = r
    d2 = sum((i - pkt[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n * n - 1))


if __name__ == '__main__':
    d = os.path.expanduser(sys.argv[1])
    maks = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    baris = []
    for tag in sorted(os.listdir(d))[:maks]:
        od = os.path.join(d, tag, 'output')
        if not os.path.isdir(od):
            continue
        t = [x for x in os.listdir(od) if x.endswith('.txt')]
        if not t:
            continue
        for gas, (bl, q, e) in blok_akhir(os.path.join(od, t[0])).items():
            if q <= 0:
                continue
            rho = spearman5(bl)
            baris.append(dict(tag=tag, gas=gas, q=q, galat=e,
                              galat_rel=100 * e / q, rho=rho,
                              hanyut='ya' if abs(rho) == 1.0 else ''))
    n = len(baris)
    hy = [b for b in baris if b['hanyut']]
    print('=' * 84)
    print('UJI KONVERGENSI BLOK — %d pengamatan dari %s' % (n, d))
    print('=' * 84)
    print('blok berurut sempurna (|rho| = 1) : %d (%.1f %%)' % (len(hy), 100 * len(hy) / n))
    print('peluang kalau murni acak          : 1,7 %')
    print()
    print('%-14s %6s %8s %8s %8s' % ('galat relatif', 'n', 'hanyut', '%', 'median |rho|'))
    tepi = [0, 5, 10, 20, 50, 100, 1e9]
    for a, b in zip(tepi, tepi[1:]):
        s = [x for x in baris if a <= x['galat_rel'] < b]
        if not s:
            continue
        h = sum(1 for x in s if x['hanyut'])
        print('%-14s %6d %8d %7.1f%% %8.2f' % (
            '%g - %g %%' % (a, b), len(s), h, 100 * h / len(s),
            st.median([abs(x['rho']) for x in s])))
    print()
    print('Kalau baris "hanyut %" tetap dekat 1,7 % di semua tingkat galat,')
    print('galat besar itu DERAU, bukan produksi yang kurang panjang.')
    print('Kalau ia naik tajam pada galat besar, produksinya memang kurang.')
