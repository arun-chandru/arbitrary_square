"""Exact continuum repair of k=5 on [643/800,129/160].

No eigenvalue approximation or sampled-parameter acceptance is used. All
acceptance tests are rational/integer and survive Python -O. The resulting
JSON is a derived certificate, not an assumed input.
"""

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

PLO, PHI = F(103993, 33102), F(104348, 33215)
ROOT_SCALE = 2**60
KAP_SPLIT = F(111, 100)
LO, HI, REF = F(643, 800), F(129, 160), F(80623, 100000)
INDEX, TAU = 5, 20


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def certify_pi():
    def atan_bounds(x):
        value = sum(((-1)**j*x**(2*j+1)/(2*j+1)
                     for j in range(12)), F(0))
        return value, value+x**25/25
    a, b = atan_bounds(F(1, 5)), atan_bounds(F(1, 239))
    require(PLO < 16*a[0]-4*b[1], 'lower Machin pi bound')
    require(16*a[1]-4*b[0] < PHI, 'upper Machin pi bound')


def sqrt_bounds(x):
    require(x >= 0, 'nonnegative radical')
    n = isqrt(x.numerator*ROOT_SCALE**2//x.denominator)
    lo, hi = F(n, ROOT_SCALE), F(n+1, ROOT_SCALE)
    require(lo*lo <= x < hi*hi, 'square-root enclosure')
    return lo, hi


def beta_bounds(kappa):
    require(1 <= kappa and kappa*kappa <= 2, 'beta domain')
    lo, hi = sqrt_bounds(2/(kappa*kappa)-1)
    return (1-hi)/2, (1-lo)/2


def phi(beta):
    require(0 < beta < 1, 'positive beta domain')
    return beta*beta/(1-beta)


def gap_counts():
    # Directly enumerate every possible positive ordered lattice pair.
    cap = isqrt(TAU)
    outer = sum(m*m+n*n < TAU
                for m in range(1, cap+1) for n in range(1, cap+1))
    inner = sum(F(m*m+n*n) <= REF*REF*TAU
                for m in range(1, cap+1) for n in range(1, cap+1))
    return outer, inner


def bbox_ratio(q):
    factor = F(1) if q >= REF else (
        q/REF*((1-KAP_SPLIT*REF)/(1-KAP_SPLIT*q))**3)
    return PLO*(1-q*q)*TAU*factor/(4*INDEX)


def angle_ratio_lower(a, b):
    require(KAP_SPLIT <= a < b <= 1/LO, 'angular box')
    require(a*LO > REF, 'reference within fixed bbox')
    bp_lower = beta_bounds(a*LO/REF)[0]
    bs_upper = beta_bounds(b)[1]
    transport = phi(bp_lower)/phi(bs_upper)
    require(0 < transport <= 1, 'transport lower bound')
    # Only q<=REF needs angular transport; q>=REF uses nesting.
    return PLO*(1-REF*REF)*TAU*transport/(4*INDEX)


def main():
    certify_pi()
    require(0 < LO < REF < HI < 1, 'target/reference order')
    require(LO*LO > F(1, 2), 'feasible angular range below sqrt(2)')
    require(KAP_SPLIT*REF < 1, 'size-growth feasibility')
    outer, inner = gap_counts()
    gap_index = outer-inner+1
    require((outer, inner, gap_index) == (11, 8, 4), 'reference gap counts')
    require(gap_index <= INDEX, 'reference indexed lower bound')

    # The area-weighted bbox factor has its minimum at an endpoint on
    # either side of REF. Its value at REF is no smaller than these.
    bbox = min(bbox_ratio(LO), bbox_ratio(HI))
    require(bbox >= 1, 'small-kappa repair')
    direct = PLO*(1-HI*HI)*TAU/(4*INDEX)
    require(direct >= 1, 'same-center nesting for q>=REF')

    stack, leaves, tests = [(KAP_SPLIT, 1/LO, 0)], [], 0
    while stack:
        a, b, depth = stack.pop()
        lower = angle_ratio_lower(a, b)
        tests += 1
        if lower >= 1:
            leaves.append((a, b, lower, depth))
        else:
            require(depth < 30, 'unresolved angular interval')
            mid = (a+b)/2
            stack.extend(((a, mid, depth+1), (mid, b, depth+1)))
    leaves.sort(key=lambda row: row[0])
    require(leaves[0][0] == KAP_SPLIT, 'left angular endpoint')
    require(leaves[-1][1] == 1/LO, 'right angular endpoint')
    for left, right in zip(leaves, leaves[1:]):
        require(left[1] == right[0], 'complete closed angular cover')
    worst = min(row[2] for row in leaves)
    summary = F(worst.numerator*10**9//worst.denominator, 10**9)
    require(summary >= 1, 'downward-rounded exact margin')

    report = {
        'status': 'PASS_EXACT_FIFTH_MODE_CONTINUUM_REPAIR',
        'scope': 'All feasible angles and all strictly interior positions; '
                 'only k=5 on the stated closed q band.',
        'q_band': [str(LO), str(HI)],
        'index': INDEX, 'reference': str(REF), 'tau': TAU,
        'strict_outer_count': outer, 'inclusive_inner_count': inner,
        'gap_index': gap_index, 'kappa_split': str(KAP_SPLIT),
        'small_kappa_ratio_lower': str(bbox),
        'direct_high_q_ratio_lower': str(direct),
        'angle_q_band': [str(LO), str(REF)],
        'angle_kappa_band': [str(KAP_SPLIT), str(1/LO)],
        'angle_leaves': len(leaves), 'angle_tests': tests,
        'maximum_depth': max(row[3] for row in leaves),
        'angle_ratio_lower_at_least': str(summary),
        'pi_lower': str(PLO), 'pi_upper': str(PHI),
        'root_scale': ROOT_SCALE,
        'dependencies': 'Uniform strict-outer/inclusive-inner gap lemma; '
                        'bbox size transport; fixed-box angle transport '
                        'L^2/U; same-center nesting.',
        'leaves': [{'kappa': [str(a), str(b)], 'ratio_lower': str(value),
                    'depth': depth} for a, b, value, depth in leaves],
    }
    out = Path(__file__).with_name('fifth_mode_upper_repair.json')
    out.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: value for key, value in report.items()
                      if key != 'leaves'}, indent=2))


if __name__ == '__main__':
    main()
