"""Independent reconstruction of the q in [12/25,7/10] angle certificate.

Reads the certificate JSON, imports neither its generator nor gap_probe,
reconstructs every lattice count and reference envelope, and checks each
closed interval and each repaired index. No files are written.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import numpy as np


def require(test, message):
    if not test:
        raise RuntimeError(message)


def floor(value):
    return value.numerator // value.denominator


def pi_bracket():
    def atan_bounds(x):
        value = sum(((-1) ** j * x ** (2*j+1) / (2*j+1)
                     for j in range(12)), F(0))
        return value, value + x ** 25 / 25
    a, b = atan_bounds(F(1, 5)), atan_bounds(F(1, 239))
    require(PLO < 16*a[0]-4*b[1], "strict lower pi bound")
    require(16*a[1]-4*b[0] < PHI, "strict upper pi bound")


CAP, DEN, S = 250000, 64000, 2**40
PLO, PHI, KAP = F(103993, 33102), F(104348, 33215), F(99, 70)


def coefficient(l, h, r):
    require(0 < l < h < 1 and 0 < r < 1, "coefficient input domain")
    require(KAP*max(h, r) < 1, "strict reference feasibility")
    def endpoint(q):
        transport = (F(1) if q >= r else
                     q*(1-KAP*r)**3 / (r*(1-KAP*q)**3))
        return (1-q*q)*transport
    return floor(S*PLO*min(endpoint(l), endpoint(h))/4)


def main():
    pi_bracket()
    require(KAP*KAP > 2, "all angles covered")
    require(S*CAP < 2**63 and DEN**2*CAP < 2**63, "vector product ranges")
    data = json.loads(Path(__file__).with_name("all_angles_certificate.json")
                      .read_text(encoding="utf-8-sig"))
    require((data["denominator"], data["scale"], data["count_cap"])
            == (DEN, S, CAP), "fixed arithmetic parameters")
    require(F(data["kappa_upper"]) == KAP, "fixed kappa upper bound")

    # Construct a sorted multiset of all squared eigenvalues. Searchsorted
    # gives the inclusive count with every ordered-pair multiplicity intact.
    modes = sorted(m*m+n*n for m in range(1, isqrt(CAP)+1)
                   for n in range(1, isqrt(CAP-m*m)+1))
    require(len(modes) == 195837, "direct lattice census")
    tau = np.arange(1, CAP+1, dtype=np.int64)
    cc = np.searchsorted(np.asarray(modes, dtype=np.int64),
                         np.arange(CAP+1, dtype=np.int64), side="right")
    for energy in (0, 1, 2, 5, 10, 40401, 100000, CAP):
        row_count = sum(isqrt(energy-m*m)
                        for m in range(1, isqrt(energy)+1))
        require(int(cc[energy]) == row_count, "independent radial count")

    def frontier(numerator, kmax):
        require(type(numerator) is int and 0 < numerator < DEN,
                "reference numerator")
        inner = numerator*numerator*tau // (DEN*DEN)
        indices = cc[tau-1]-cc[inner]+1
        require(bool(np.all(indices >= 1)), "positive indexed gaps")
        envelope = np.zeros(kmax+1, dtype=np.int64)
        used = indices <= kmax
        np.maximum.at(envelope, indices[used], tau[used])
        np.maximum.accumulate(envelope, out=envelope)
        return envelope[1:]

    records = data["records"]
    require(len(records) == 315, "declared interval count")
    previous = F(12, 25)
    pairs = repairs = maximum = 0
    for rec in records:
        l, h = map(F, rec["interval"])
        require(l == previous and l < h <= F(7, 10), "closed interval coverage")
        previous = h
        kmax = floor(PHI*(1-l*l)*201**2/(4*l*l))
        require(rec["kmax"] == kmax and 1 <= kmax < CAP,
                "conservative finite index cutoff")
        ks = np.arange(1, kmax+1, dtype=np.int64)
        passed = ks <= 2  # Faber--Krahn and Krahn--Szego, analytic dependencies.
        require(len(rec["references"]) == 2, "two default references")
        for ref in rec["references"]:
            numerator = ref["numerator"]
            mu = coefficient(l, h, F(numerator, DEN))
            require(ref["coefficient"] == mu and 0 <= mu < S,
                    "exact conservative default coefficient")
            passed |= mu*frontier(numerator, kmax) >= S*ks

        gram = l >= F(3, 5) and PLO*(1-h*h)*F(1117, 109) >= 20
        require(rec["gram_for_k5"] is gram, "exact fifth-mode applicability")
        if gram:
            passed[4] = True
        missing = ks[~passed].tolist()
        require(rec["initial_missing"] == missing, "reconstructed missing list")
        repaired = []
        for fix in rec["repairs"]:
            k, numerator, energy = fix["k"], fix["reference"], fix["tau"]
            require(all(type(v) is int for v in (k, numerator, energy)),
                    "integer repair data")
            require(k in missing and 1 <= energy <= CAP, "repair target and energy")
            mu = coefficient(l, h, F(numerator, DEN))
            inner = numerator*numerator*energy // (DEN*DEN)
            gap_index = int(cc[energy-1]-cc[inner])+1
            require(fix["coefficient"] == mu, "exact repair coefficient")
            require(1 <= gap_index == fix["gap_index"] <= k,
                    "strict outer / inclusive inner repair gap")
            require(mu*energy >= S*k, "repaired Polya inequality")
            repaired.append(k)
        require(sorted(repaired) == missing, "all missing indices repaired exactly once")
        pairs += kmax
        repairs += len(repaired)
        maximum = max(maximum, kmax)
    require(previous == F(7, 10), "final endpoint")
    require((pairs, repairs, maximum) == (17944637, 28, 105989),
            "independently reproduced census")
    require(data["pairs"] == pairs and data["repairs"] == repairs
            and data["max_index"] == maximum and data["unresolved"] == [],
            "saved summary agrees with reconstructed proof")
    print(json.dumps({"status": "PASS_INDEPENDENT_REPLAY", "intervals": len(records),
                      "pairs": pairs, "repairs": repairs, "max_index": maximum,
                      "scope": "finite indexed part on [12/25,7/10]",
                      "analytic_dependencies": ["rotated bounding-box transport",
                       "strict/inclusive endpoint gap", "q<=4/5 disk corridor tail",
                       "Faber--Krahn and Krahn--Szego", "uniform fifth-mode Gram"]},
                     indent=2))


if __name__ == "__main__":
    main()
