"""Exact continuum repairs for the ninth and eleventh rotated-hole modes.

All acceptance comparisons use Fraction and integer square roots. No floating
arithmetic, eigenvalue approximation, or sampled-parameter claim is used.
Run with Python -B -O. Only the derived JSON beside this file is written.
"""

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

PLO, PHI = F(103993, 33102), F(104348, 33215)
ROOT_SCALE = 2**60
KAP_SPLIT = F(111, 100)
BANDS = (
    (9, F(299, 400), F(3, 4), F(3, 4), 32),
    (11, F(3153, 4000), F(3163, 4000), F(79057, 100000), 40),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def certify_pi():
    def atan_bounds(x):
        val = sum(((-1)**j * x**(2*j+1) / (2*j+1)
                   for j in range(12)), F(0))
        following = x**25 / 25
        return val, val + following
    a, b = atan_bounds(F(1, 5)), atan_bounds(F(1, 239))
    require(PLO < 16*a[0]-4*b[1], 'strict lower pi bound')
    require(16*a[1]-4*b[0] < PHI, 'strict upper pi bound')


def sqrt_bounds(x):
    require(x >= 0, 'nonnegative radicand')
    n = isqrt(x.numerator * ROOT_SCALE**2 // x.denominator)
    lo, hi = F(n, ROOT_SCALE), F(n+1, ROOT_SCALE)
    require(lo*lo <= x < hi*hi, 'verified radical enclosure')
    return lo, hi


def beta_bounds(kappa):
    require(1 <= kappa and kappa*kappa <= 2, 'beta argument range')
    lo, hi = sqrt_bounds(2/(kappa*kappa)-1)
    return (1-hi)/2, (1-lo)/2


def f(beta):
    require(0 < beta < 1, 'positive beta in f')
    return beta*beta/(1-beta)


def gap_index(r, tau):
    # Independent direct enumeration of positive ordered lattice pairs.
    cap = isqrt(tau)
    outer = sum(m*m+n*n < tau
                for m in range(1, cap+1) for n in range(1, cap+1))
    inner = sum(F(m*m+n*n) <= r*r*tau
                for m in range(1, cap+1) for n in range(1, cap+1))
    return outer-inner+1, outer, inner


def box_transport(q, r):
    if r <= q:
        return F(1)
    return q/r * ((1-KAP_SPLIT*r)/(1-KAP_SPLIT*q))**3


def angle_ratio_lower(l, h, r, tau, k, a, b):
    require(KAP_SPLIT <= a < b <= 1/l, 'kappa box range')
    require(l <= h <= r and a*l > r, 'angle-reduction reference range')
    bp_lo = beta_bounds(a*l/r)[0]
    bs_hi = beta_bounds(b)[1]
    transport = f(bp_lo)/f(bs_hi)
    require(0 < transport <= 1, 'angle transport lower bound')
    return PLO*(1-h*h)*tau*transport/(4*k)


def main():
    certify_pi()
    results = []
    for k, l, h, r, tau in BANDS:
        require(0 < l <= h < 1 and l*l > F(1, 2), 'target band')
        require(KAP_SPLIT*r < 1, 'uniform size-growth feasibility')
        idx, outer, inner = gap_index(r, tau)
        require(idx <= k, 'reference indexed gap')

        # kappa <= KAP_SPLIT: monotonicity in kappa and the proved
        # endpoint minimum in q bound the entire real rectangle.
        bbox = min(PLO*(1-q*q)*box_transport(q, r)*tau/(4*k)
                   for q in (l, h))
        require(bbox >= 1, 'small-kappa continuum repair')

        # q >= r requires only a same-center smaller reference hole.
        direct = None
        if h >= r:
            direct = PLO*(1-h*h)*tau/(4*k)
            require(direct >= 1, 'direct high-q repair')

        ah = min(h, r)
        leaves = []
        stack = [(KAP_SPLIT, 1/l, 0)]
        while stack:
            a, b, depth = stack.pop()
            lower = angle_ratio_lower(l, ah, r, tau, k, a, b)
            if lower >= 1:
                leaves.append((a, b, lower, depth))
            else:
                require(depth < 30, 'unresolved angle interval')
                mid = (a+b)/2
                stack.extend(((a, mid, depth+1), (mid, b, depth+1)))

        leaves.sort(key=lambda row: row[0])
        require(leaves[0][0] == KAP_SPLIT, 'left kappa endpoint')
        require(leaves[-1][1] == 1/l, 'right kappa endpoint')
        for left, right in zip(leaves, leaves[1:]):
            require(left[1] == right[0], 'no angular coverage gap')
        worst = min(row[2] for row in leaves)
        # The exact proof keeps every Fraction; the summary rounds down.
        summary = F(worst.numerator*10**9//worst.denominator, 10**9)
        require(summary >= 1, 'summary direction and margin')
        results.append({
            'k': k, 'q_band': [str(l), str(h)], 'reference': str(r),
            'tau': tau, 'gap_index': idx,
            'strict_outer_count': outer, 'inclusive_inner_count': inner,
            'kappa_split': str(KAP_SPLIT),
            'small_kappa_ratio_lower': str(bbox),
            'direct_high_q_ratio_lower': str(direct) if direct else None,
            'angle_q_band': [str(l), str(ah)],
            'angle_kappa_band': [str(KAP_SPLIT), str(1/l)],
            'angle_intervals': len(leaves),
            'maximum_depth': max(row[3] for row in leaves),
            'angle_ratio_lower_at_least': str(summary),
            'leaves': [{'kappa': [str(a), str(b)],
                        'ratio_lower': str(value), 'depth': depth}
                       for a, b, value, depth in leaves],
        })

    report = {
        'status': 'PASS_EXACT_CONTINUUM_REPAIRS',
        'scope': 'All feasible angles and all strictly interior placements '
                 'for the two stated q/index bands.',
        'dependencies': 'Strict-outer/inclusive-inner gap lemma; '
                        'rotated bounding-box size transport; '
                        'fixed-box rotation transport L^2/U; same-center '
                        'domain monotonicity.',
        'pi_lower': str(PLO), 'pi_upper': str(PHI),
        'root_scale': ROOT_SCALE, 'records': results,
    }
    out = Path(__file__).with_name('low_mode_rotation_repairs.json')
    out.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({
        'status': report['status'],
        'records': [{key: value for key, value in row.items() if key != 'leaves'}
                    for row in results],
    }, indent=2))


if __name__ == '__main__':
    main()
