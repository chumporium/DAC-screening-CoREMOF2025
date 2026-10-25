"""Classify K_H values by whether they are converged, not merely by a small error bar.

A small error bar does not imply a correct value. In heavy-tailed systems all
blocks can agree on the same wrong value as long as rare events have not yet
been sampled. The distinguishing criteria are whether the value still drifts
when the run is extended, and whether the mean is set by a handful of events.

Three classes:

  KONVERGEN        (converged) blocks are consistent and the value does not
                   shift when the number of cycles is increased. Valid value.

  KONVERGEN_LAMBAT (slowly converging) heavy-tailed at the short run length but
                   levels off at the long run length. The short-run value is
                   wrong; the long-run value is valid.

  TAK_KONVERGEN    (not converged) a single block dominates regardless of run
                   length. The cumulative mean decays as 1/k, the signature of a
                   single event divided among more and more blocks rather than
                   convergence to any value. Widom insertion cannot measure
                   these materials; the value is an artefact.

The 10x dominance threshold is deliberately loose: passing this test does not
mean the estimate is good, only that it is not a single-event artefact.
"""
import os, re, sys, csv, io, glob, math, statistics as st

T95 = 2.7764
DOMINAN = 10.0          # max/median in the long run; above this = single event
GESER = 1.5             # acceptable shift of the value between the two run lengths

PB = re.compile(r'Henry coefficient based on Rosenbluth weight:(.*?)Average Henry', re.S)
PBLK = re.compile(r'Block\[\s*\d+\]\s+(-?\d[\d.eE+-]*)')
PAT = re.compile(r'Average Henry coefficient:\s+([\d.eE+-]+)\s+\+/-\s+([\d.eE+-]+)\s+\[mol/kg/Pa\]')


def baca(d):
    """-> (KH, sigma_pct, blocks[5]) or None."""
    f = glob.glob(os.path.join(d, 'output', '*.txt'))
    if not f:
        return None
    t = open(f[0], encoding='utf-8', errors='replace').read()
    m = PAT.search(t)
    if not m:
        return None
    kh, ci = float(m.group(1)), float(m.group(2))
    if kh <= 0:
        return None
    mb = PB.search(t)
    b = []
    if mb:
        for y in PBLK.findall(mb.group(1)):
            try:
                b.append(float(y))
            except ValueError:
                pass
    return kh, 100 * ci / T95 / kh, (b if len(b) == 5 else None)


def luruh_1k(b):
    """Check whether the cumulative mean decays as 1/k.

    If a single block contributes everything and the others are zero, the
    cumulative mean after k blocks is S/k, so kum[k]*k is nearly constant.
    That is not convergence, only one event divided among more blocks.
    """
    if not b or min(b) < 0:
        return False
    kum = [st.mean(b[:k + 1]) for k in range(5)]
    if kum[0] <= 0:
        return False
    prod = [kum[k] * (k + 1) for k in range(5)]
    # if the product is nearly constant (spread < 30 %), 1/k pattern
    return (max(prod) / min(prod)) < 1.3 if min(prod) > 0 else False


def golongkan(pendek, panjang):
    """pendek/panjang: output of baca() for the short and long run."""
    if panjang is None:
        return 'BELUM_DIUJI', None
    khp, sgp, bp = panjang
    mm = (max(bp) / st.median(bp)) if bp and st.median(bp) > 0 else float('nan')
    if bp and (mm > DOMINAN or luruh_1k(bp)):
        return 'TAK_KONVERGEN', mm
    if pendek is None:
        return 'KONVERGEN', mm
    kh0 = pendek[0]
    nis = khp / kh0 if kh0 > 0 else float('inf')
    if nis > GESER or nis < 1 / GESER:
        return 'KONVERGEN_LAMBAT', mm
    return 'KONVERGEN', mm


if __name__ == '__main__':
    K20 = {}
    p = os.path.expanduser('~/ulangco2_hasil.csv')
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            K20[r['nama']] = (float(r['KH']) if 'KH' in r else float(r['KH_baru']),
                              float(r['sigma_pct']), None)

    baris = []
    # long runs: ~/lonjak (500k cycles) and ~/konv (200k, 2M cycles)
    for base, sik in [(os.path.expanduser('~/lonjak'), 500000)]:
        if not os.path.isdir(base):
            continue
        for n in sorted(os.listdir(base)):
            r = baca(os.path.join(base, n))
            if r is None:
                continue
            gol, mm = golongkan(K20.get(n), r)
            baris.append(dict(nama=n, siklus=sik, KH='%.6e' % r[0],
                              sigma_pct='%.2f' % r[1],
                              maks_per_median=('%.2f' % mm) if mm == mm else '',
                              golongan=gol))
    base = os.path.expanduser('~/konv')
    if os.path.isdir(base):
        for d in sorted(os.listdir(base)):
            if '__s' not in d:
                continue
            n, sik = d.split('__s')
            r = baca(os.path.join(base, d))
            if r is None:
                continue
            gol, mm = golongkan(K20.get(n), r)
            baris.append(dict(nama=n, siklus=int(sik), KH='%.6e' % r[0],
                              sigma_pct='%.2f' % r[1],
                              maks_per_median=('%.2f' % mm) if mm == mm else '',
                              golongan=gol))
    if not baris:
        raise SystemExit('belum ada lari panjang yang selesai')
    out = os.path.expanduser('~/golongan_konvergensi.csv')
    with io.open(out, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(baris[0].keys()))
        w.writeheader()
        for b in baris:
            w.writerow(b)
    import collections
    c = collections.Counter(b['golongan'] for b in baris)
    print('%d lari panjang dinilai' % len(baris))
    for k, v in c.most_common():
        print('  %-18s %d' % (k, v))
    print('->', out)
