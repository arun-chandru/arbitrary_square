"""Independent exact replay of every q in [7/10,4/5] finite-index box.

Imports no generator, gap_probe, or other replay. Reconstructs the lattice
multiset and checks each target index by a suffix minimum over admissible
future energies, independent of the generator's cumulative-max frontier.
The only reported exceptions must lie inside the two proved low-mode bands.
No files are written. Compatible with Python -B -O.
"""

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import numpy as np

CAP, SCALE = 250000, 2**40
PLO, PHI = F(103993, 33102), F(104348, 33215)
REPAIR_BANDS = {9: (F(299, 400), F(3, 4)),
                11: (F(3153, 4000), F(3163, 4000))}


def require(test, message):
    if not test:
        raise RuntimeError(message)


def floor(value):
    return value.numerator//value.denominator


def certify_pi():
    def atan(x):
        value = sum(((-1)**j*x**(2*j+1)/(2*j+1)
                     for j in range(12)), F(0))
        return value, value+x**25/25
    a, b = atan(F(1, 5)), atan(F(1, 239))
    require(PLO < 16*a[0]-4*b[1] < 16*a[1]-4*b[0] < PHI,
            'exact rational Machin bracket')


def main():
    certify_pi()
    require(SCALE*CAP < 2**63 and 64000**2*CAP < 2**63,
            'integer array product ranges')
    # A sorted list keeps every ordered pair, including multiplicities.
    modes = np.asarray(sorted(m*m+n*n
                       for m in range(1, isqrt(CAP)+1)
                       for n in range(1, isqrt(CAP-m*m)+1)), dtype=np.int64)
    require(len(modes) == 195837, 'complete disk census')
    energies = np.arange(CAP+1, dtype=np.int64)
    count = np.searchsorted(modes, energies, side='right')
    for e in (0, 1, 2, 5, 18, 25, 32, 40, 40401, 100000, CAP):
        direct = sum(isqrt(e-m*m) for m in range(1, isqrt(e)+1))
        require(int(count[e]) == direct, 'independent row census at '+str(e))
    tau = energies[1:]

    data = json.loads(Path(__file__).with_name('certificate_midband.json')
                      .read_text(encoding='utf-8-sig'))
    records = data['records']
    require(len(records) == 267, 'fixed interval count')
    previous = F(7, 10)
    pairs = maximum = 0
    unmatched = []
    repair_intervals = {9: [], 11: []}
    for record in records:
        l, h = map(F, record['interval'])
        require(previous == l < h <= F(4, 5), 'closed size interval coverage')
        previous = h
        require(l.denominator <= 64000 and l.numerator < l.denominator,
                'reference arithmetic range')
        kmax = floor(PHI*(1-l*l)*201**2/(4*l*l))
        require(record['kmax'] == kmax and 2 <= kmax < CAP,
                'complete finite index cutoff')
        mu = floor(SCALE*PLO*(1-h*h)/4)
        require(record['mu'] == mu and 0 < mu < SCALE,
                'downward area/pi coefficient')

        inner = (l.numerator*l.numerator*tau)//(l.denominator*l.denominator)
        indices = count[tau-1]-count[inner]+1
        require(bool(np.all(indices >= 1)), 'positive reference gap indices')
        # For target k, tau must be at least ceil(SCALE*k/mu). If some
        # such tau has gap_index<=k, min--max and mu*tau>=SCALE*k prove it.
        best_index_from = np.minimum.accumulate(indices[::-1])[::-1]
        ks = np.arange(1, kmax+1, dtype=np.int64)
        required_tau = (SCALE*ks+mu-1)//mu
        require(bool(np.all((required_tau >= 1) & (required_tau <= CAP))),
                'all candidate thresholds inside exact census')
        passed = (ks <= 2) | (best_index_from[required_tau-1] <= ks)
        missing = ks[~passed].tolist()
        require(missing == record['remaining'], 'independent missing-index list')
        if missing:
            unmatched.append(record)
        for k in missing:
            require(k in REPAIR_BANDS, 'only proved repair indices remain')
            a, b = REPAIR_BANDS[k]
            require(a <= l < h <= b, 'whole failed interval inside repair band')
            repair_intervals[k].append((l, h))
        pairs += kmax
        maximum = max(maximum, kmax)

    require(previous == F(4, 5), 'final size endpoint')
    require(pairs == 6343627, 'every finite index/interval pair replayed')
    require(data['intervals'] == len(records) and data['pairs'] == pairs,
            'saved census agrees')
    require(data['remaining_indices'] == [9, 11]
            and data['remaining'] == unmatched, 'saved residual data agrees')
    for k, spans in repair_intervals.items():
        require(spans[0][0] == REPAIR_BANDS[k][0]
                and spans[-1][1] == REPAIR_BANDS[k][1], 'exact repair band endpoints')
        for x, y in zip(spans, spans[1:]):
            require(x[1] == y[0], 'contiguous repaired residual band')

    print(json.dumps({
        'status': 'PASS_INDEPENDENT_MIDBAND_REPLAY',
        'intervals': len(records), 'pairs': pairs, 'maximum_index': maximum,
        'unresolved_after_proved_rotation_repairs': [],
        'repaired_interval_index_pairs': sum(len(v) for v in repair_intervals.values()),
        'repair_bands': {str(k): [str(a), str(b)]
                         for k, (a, b) in REPAIR_BANDS.items()},
        'independent_method': 'sorted ordered lattice modes; inclusive search; '
                              'suffix-minimum gap-index verification at every k',
        'analytic_tail': 'q<=4/5 and sqrt(E)/pi>=201/q; '
                         'the existing exact disk-deficit corridor',
    }, indent=2))


if __name__ == '__main__':
    main()
