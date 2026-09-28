"""Exact rational checks for the explicit circle estimate in manuscript.tex.

This verifies only the displayed scalar comparisons, not the Bessel,
Stieltjes, or Poisson arguments. No lattice scan is performed.
"""

from fractions import Fraction as F


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


pi_lo = F(333, 106)
pi_hi = F(355, 113)
require(F(32, 13) / pi_lo < F(64, 81), "global Bessel factor")
require(F(9, 25) * 2 < F(64, 81), "small-argument Bessel interval")
factor = F(100001, 100000)
require(factor**2 * (1 - F(1, 48 * 63**2)) > 1,
        "large disk Fourier factor")

sqrt2_hi = F(283, 200)
require(sqrt2_hi**2 > 2, "sqrt2 upper bound")
require(8 * F(119, 200)**4 > 1, "2^(-3/4) upper bound")
require(3 - F(23, 4) / sqrt2_hi < 0, "C2 pi coefficient")
c2_upper = (4 + 4 * F(119, 200) + sqrt2_hi + 3 * F(157, 50)
            - (23 * F(157, 50) + 26) / (4 * sqrt2_hi))
require(c2_upper < 0, "C2 sign")

require(3200 > 81 * pi_hi**3, "c > 1/10")
require(F(32, 15) > 2, "K >= 2 at radius64")
c13_hi = F(4833, 10000)
require(2 * pi_lo**3 * c13_hi**6 > F(8, 9)**2,
        "c^(1/3) upper bound")
sqrt83_hi = F(1633, 1000)
require(sqrt83_hi**2 > F(8, 3), "sqrt(8/3) upper bound")
total = (F(3, 4) * pi_hi + F(9, 64 * 256) * pi_hi
         + 6 * factor * c13_hi * sqrt83_hi * (1 + F(3, 4096)))
require(total < F(71, 10), "final 7.1 bound")

for eta, v in [(F(2, 25), 23), (F(1, 14), 25),
               (F(1, 18), 32), (F(1, 19), 34),
               (F(1, 20), 36), (F(1, 24), 43)]:
    require(eta >= F(3, 4 * v**3) + F(71, 40 * v),
            f"corridor eta={eta}, R={v**3}")

print("PASS: all exact rational circle-tail constants and table onsets")
print("C2 rational upper bound:", c2_upper)
print("Final rational upper bound:", total)
