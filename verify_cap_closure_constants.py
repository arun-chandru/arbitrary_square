"""Exact scalar audit of the q>=49/50 conformal-cap proof.

This verifies rational budget arithmetic, not the spectral comparisons.
Those comparisons are proved in the conformal-cap and sector sections
of the accompanying standalone manuscript.
"""

from fractions import Fraction as Q
from math import factorial


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)
    print("PASS", label)


def main():
    delta = Q(1, 50)
    area = Q(2, 5) * Q(51, 50)**2 * Q(53, 50)**3 + Q(2, 5) * Q(21, 20)**3
    linear = (Q(51, 100) + 1) * 4 * Q(51, 50) + 4 * Q(1001, 1000)
    channels = 2 * Q(51, 50) * Q(53, 50) + Q(21, 20) * Q(1001, 1000)
    require(area < Q(24, 25), "true-area inflation < (24/25) delta^5")
    require(linear < Q(1017, 100), "linear corner cost < (1017/100) delta")
    require(channels < Q(161, 50), "channel coefficient <161/50")
    require(Q(22, 7) / 4 * Q(24, 25) * 216 < 163, "area-tail coefficient <163")
    require(Q(22, 7) * Q(51, 50) * 216 < 700, "top logarithm argument")
    require(Q(22, 7) * Q(1001, 1000) * 216 < 700, "right logarithm argument")
    exp3_lower = sum((Q(3)**j / factorial(j) for j in range(10)), Q(0))
    require(exp3_lower > 20, "finite Taylor lower bound exp(3)>20")
    require(20**5 > 700 * 50**2, "log(1750000)<15")
    corner = Q(1017, 100)*delta + 163*delta**2 + Q(161, 50)*delta**2*(Q(1, 2)+Q(15, 6))
    require(corner == Q(17029, 62500), "corner total =0.272464")
    require(corner < Q(7, 25), "corner total <7/25")
    require(Q(7, 20)*Q(97, 100)-Q(51, 2500) > Q(31, 100), "active sector budget >=31/100")
    require(Q(97, 100)**2 - delta**2/2 > Q(47, 50), "sector quadratic coefficient")
    require(Q(97, 100)*Q(99, 100) > Q(47, 50), "sector mixed coefficient")
    require(Q(47*3, 400) > Q(7, 20), "inactive endpoint budget >7/20 using pi>3")
    require(2*Q(31, 100) > Q(7, 20), "active endpoint budget >7/20")
    print("ALL EXACT SCALAR CHECKS PASS; operator proof remains a separate dependency.")


if __name__ == "__main__":
    main()
