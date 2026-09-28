"""Exact scalar certificate for a uniform fifth-eigenvalue hole-mass bound.

No approximate perforated-domain eigenvalues are computed. All acceptance
decisions use fractions and unconditional checks; python -O is supported.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def atan_bounds(x, n=12):
    terms = [(-1) ** j * x ** (2 * j + 1) / (2 * j + 1)
             for j in range(n)]
    partial = sum(terms, F(0))
    next_term = (-1) ** n * x ** (2 * n + 1) / (2 * n + 1)
    return min(partial, partial + next_term), max(partial, partial + next_term)


a5 = atan_bounds(F(1, 5))
a239 = atan_bounds(F(1, 239))
machin_lo = 16 * a5[0] - 4 * a239[1]
machin_hi = 16 * a5[1] - 4 * a239[0]
PLO, PHI = F(333, 106), F(355, 113)
require(PLO < machin_lo < machin_hi < PHI, "Machin pi bracket")
PMID, PRAD = (PLO + PHI) / 2, (PHI - PLO) / 2


@lru_cache(maxsize=None)
def trig_pi(a, cosine=False):
    """Enclose sin(pi*a) or cos(pi*a) using Taylor + Lipschitz pi error."""
    x = PMID * a
    require(abs(x) <= 3, "Taylor argument exceeds declared compact range")
    if cosine:
        val = sum(((-1) ** j * x ** (2 * j) / factorial(2 * j)
                   for j in range(13)), F(0))
        error = abs(x) ** 25 / factorial(25)
    else:
        val = sum(((-1) ** j * x ** (2 * j + 1) / factorial(2 * j + 1)
                   for j in range(13)), F(0))
        error = abs(x) ** 26 / factorial(26)
    error += PRAD * abs(a)
    return val - error, val + error


def lower_c(l, u):
    # A is the least sin(pi*x)^2 interval mass allowed by the center bound.
    m = (l + u) / 2
    sin_lo, sin_hi = trig_pi(m)
    sin_lo -= PHI * (u - l) / 2
    sin_hi += PHI * (u - l) / 2
    cos_lo = trig_pi(1 - F(9, 25) / u, True)[0]
    require(sin_lo > 0 and cos_lo > 0, "positive factors in A")
    a_lo = l / 2 + sin_lo * cos_lo / (2 * PHI)

    # D(t)=t/2+sin(3*pi*t)/(6*pi) is increasing. On this range,
    # sin(3*pi*l)=-sin(pi*(2-3*l)), whose latter sine is positive.
    d_lo = l / 2 - trig_pi(2 - 3 * l)[1] / (6 * PLO)
    require(d_lo > 0, "positive D")

    # |sin(2*pi*t)|=sin(pi*|2*t-1|), increasing in |2*t-1| here.
    z = max(abs(2 * l - 1), abs(2 * u - 1))
    b_hi = sin_hi / (2 * PLO) + trig_pi(z)[1] / (4 * PLO)
    require(b_hi > 0, "positive B")
    return 4 * (a_lo * d_lo - b_hi * b_hi)


def main():
    # [53/125,3/5] contains [3/(5*sqrt(2)),3/5].
    require(F(53, 125) ** 2 < F(9, 50), "irrational left endpoint coverage")
    least = None
    worst = None
    for n in range(424, 600):
        l, u = F(n, 1000), F(n + 1, 1000)
        value = lower_c(l, u)
        require(value >= F(9, 100), f"Unresolved interval [{l},{u}]")
        if least is None or value < least:
            least, worst = value, (l, u)
    require(least is not None, "nonempty coverage")
    # A short rational summary is independently weaker than the exact minimum.
    summary_lower = F(least.numerator * 10**6 // least.denominator, 10**6)
    require(summary_lower >= F(9, 100), "summary rounding direction")
    print(json.dumps({
        "status": "PASS",
        "intervals": 176,
        "covered": ["53/125", "3/5"],
        "required": "9/100",
        "certified_lower_at_least": str(summary_lower),
        "worst_interval": [str(v) for v in worst],
        "eigenvalue_consequence": "lambda_5 >= (1117/109)*pi^2 for q>=3/5",
        "polya_consequence": "q^2 <= 1 - 2180/(1117*pi)",
        "arithmetic": "fractions; Taylor remainders; rational Machin pi bounds",
    }, indent=2))


if __name__ == "__main__":
    main()
