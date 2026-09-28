"""Exact scalar checks for the manuscript's q>=23/25, X>=8 cap bound.

No spectral scan is performed; the operator comparisons are proved in
the accompanying standalone manuscript.
"""

from fractions import Fraction as F


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


delta = F(2, 25)
pi_hi = F(22, 7)
tail = F(5832, 125)
require(1 - F(1, 200) + F(1, 10) - F(1, 6000) > F(25, 23),
        "theta < 1/10")
require(1 - F(1, 20) - F(1, 600) > F(237, 250),
        "kappa Taylor coefficient")
require(F(23, 25) * F(237, 250) * F(7, 6) > 1,
        "theta <= 7delta/6")
require(F(7, 6) * F(200, 199) < F(6, 5), "tangent bound")
require(1 + F(49, 72) * delta < F(53, 50), "w bound")
require(1 - F(53, 50) * delta - F(3, 5) * delta**2 > F(91, 100),
        "sector length")
require(F(400, 399) < F(1003, 1000), "sqrt gamma")

z = F(1, 100)
require(F(1, 100) - F(77, 400) * z > 0,
        "gamma quadratic bound over interval")
require((1 + F(3, 200) * z)**2 * (1 - z / 36) > 1,
        "beta bound at endpoint")
require(F(51, 100) + F(3, 200) + F(153, 20000) * z < F(3, 5),
        "gamma/beta inflation")
require(pi_hi * F(53, 50) * tail < 160, "all logarithm arguments")

area_hi = F(2, 5) * (F(53, 50)**2 * F(6, 5)**3 + F(7, 6)**3)
require(area_hi < F(71, 50), "high angle area inflation")
require(pi_hi * tail * F(71, 50) / 4 < F(521, 10),
        "high angle area count")
channels_hi = 2 * F(53, 50) * F(6, 5) + F(7, 6) * F(1003, 1000)
require(channels_hi < F(93, 25), "high angle channels")
budget_hi = 4 * (F(7, 20) * F(91, 100) - F(53, 50) * delta)
require(F(3, 10) * F(91, 100)**2 > budget_hi / 4,
        "high angle inactive sector budget")
error_hi = (F(141139, 25000) * delta + F(521, 10) * delta**2
            + F(93, 25) * delta**2 * F(7, 3) + 3 * delta / 8)
require(error_hi < budget_hi, "high angle total")

require(1 / (1 - F(1, 2 * 125**2)) < F(1001, 1000),
        "low angle tangent multiplier")
require(1 + delta / 200 < F(1001, 1000), "low angle w")
require(F(23, 25) - (F(1, 200) + F(1001, 20000)) * delta**2
        > F(919, 1000), "low angle sector length")
area_lo = F(2, 5) * (F(1001, 1000)**2 * F(1001, 10000)**3 + F(1, 1000))
require(area_lo < F(803, 1000000), "low angle area inflation")
require(pi_hi * tail * F(803, 1000000) / 4 < F(3, 100),
        "low angle area count")
channels_lo = 2 * F(1001, 1000) * F(1001, 10000) + F(1, 10) * F(1001, 1000)
require(channels_lo < F(301, 1000), "low angle channels")
budget_lo = 2 * (F(7, 20) * F(919, 1000) - F(1001, 1000) * delta)
error_lo = (F(2661, 500) * F(1001, 1000) * delta + F(3, 100) * delta**2
            + F(301, 1000) * delta**2 * F(7, 3) + 3 * delta / 8)
require(error_lo < budget_lo, "low angle total")

print("PASS: all displayed q>=23/25, X>=8 scalar comparisons")
print("High-angle error and deficit:", error_hi, budget_hi)
print("Low-angle error and deficit:", error_lo, budget_lo)
