"""Exact compact certificate for D(y)>=7/20 on 1<=y<=5/2.

The manuscript supplies the analytic tail for y>=5/2. This verifier uses
75 closed intervals, rational Machin bounds, and integer square roots.
Run with Python -B -O; no data files or floating-point arithmetic are used.
"""

from fractions import Fraction as F
from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    def atan_bounds(x):
        value = sum(((-1)**j*x**(2*j+1)/(2*j+1)
                     for j in range(12)), F(0))
        return value, value+x**25/25
    a, b = atan_bounds(F(1, 5)), atan_bounds(F(1, 239))
    pi_lower = F(333, 106)
    require(pi_lower < 16*a[0]-4*b[1], 'Machin lower bound')
    scale = 2**50

    def sqrt_upper(x):
        require(x >= 0, 'nonnegative radicand')
        n = isqrt(x.numerator*scale**2//x.denominator)
        lo, hi = F(n, scale), F(n+1, scale)
        require(lo*lo <= x < hi*hi, 'square-root enclosure')
        return hi

    worst = None
    previous = F(1)
    for j in range(75):
        left, right = F(1)+F(j, 50), F(1)+F(j+1, 50)
        require(left == previous, 'closed interval coverage')
        previous = right
        upper_radicals = sum((sqrt_upper(1-F(m*m)/(right*right))
                              for m in range(1, right.numerator//right.denominator+1)),
                             F(0))
        lower = pi_lower*left/4-upper_radicals
        require(lower >= F(7, 20), 'compact radial deficit')
        worst = lower if worst is None else min(worst, lower)
    require(previous == F(5, 2), 'final endpoint')
    # The elementary tail is D>=1/2-1/sqrt(18y). At y>=5/2,
    # 1/(18y)<=1/45<9/400=(3/20)^2.
    require(F(1, 45) < F(9, 400), 'analytic tail threshold')
    rounded = F(worst.numerator*10**8//worst.denominator, 10**8)
    require(rounded >= F(7, 20), 'downward summary')
    print('PASS_RADIAL_SCALAR: 75 closed intervals; exact minimum lower '
          'bound at least '+str(rounded)+'; analytic tail from 5/2.')


if __name__ == '__main__':
    main()
