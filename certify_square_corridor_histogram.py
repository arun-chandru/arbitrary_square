"""Exact all-real square-count corridors from rational-radius bins.

This is NOT an eigenvalue discretization. Every generated lattice pair is
placed using a checked integer inequality. A rational polynomial envelope
certifies all radii between adjacent bins. Analytic extension beyond R is a
separate dependency. Memory is O(h*R), not O(R**2); generation is single-process.
"""

import argparse
from fractions import Fraction as Q
import json
from math import isqrt
from pathlib import Path
import time

import numpy as np


def atan_interval(x, n):
    s = sum(((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n)), Q(0))
    nxt = (-1)**n*x**(2*n+1)/Q(2*n+1)
    return min(s, s+nxt), max(s, s+nxt)


def pi_intervals():
    a, b = atan_interval(Q(1, 5), 24)
    c, d = atan_interval(Q(1, 239), 8)
    lo, hi = 16*a-4*d, 16*b-4*c
    low, high = Q(103993, 33102), Q(104348, 33215)
    if not low < lo < hi < high:
        raise ArithmeticError("Machin pi enclosure")
    return low, high


def direct_count(j, h):
    # h^2(m^2+n^2)<=j^2 iff m^2+n^2<=floor(j^2/h^2).
    bound = j*j//(h*h)
    return sum(isqrt(bound-m*m) for m in range(1, isqrt(bound)+1))


def build_histogram(R, h):
    size = R*h
    if size*size*104348 >= 2**63-1:
        raise ValueError("requested radius exceeds audited int64 products")
    hist = np.zeros(size+1, dtype=np.int64)
    pairs = 0
    start = time.monotonic()
    next_notice = start+20
    for m in range(1, isqrt(R*R//2)+1):
        upper = isqrt(R*R-m*m)
        n = np.arange(m, upper+1, dtype=np.int64)
        squared = (n*n+m*m)*(h*h)
        idx = np.ceil(np.sqrt(squared)).astype(np.int64)
        idx += (idx*idx < squared)
        idx -= ((idx-1)*(idx-1) >= squared)
        if np.any(idx*idx < squared) or np.any((idx-1)*(idx-1) >= squared):
            raise ArithmeticError("integer bin certificate failed")
        np.add.at(hist, idx, 2)
        hist[int(idx[0])] -= 1  # diagonal m=n counted once
        pairs += len(n)
        now = time.monotonic()
        if now >= next_notice:
            print(json.dumps({"phase":"lattice histogram", "row":m,
                              "last_row":isqrt(R*R//2), "pairs":pairs,
                              "seconds":round(now-start, 2)}), flush=True)
            next_notice = now+20
    np.cumsum(hist, out=hist)
    # Different O(R) implementation checks deterministic dispersed radii.
    check_indices = sorted({1, h, size, *[1+(size-1)*i//31 for i in range(32)]})
    for j in check_indices:
        if hist[j] != direct_count(j,h):
            raise ArithmeticError(("independent direct count",j,int(hist[j]),direct_count(j,h)))
    return hist, pairs, time.monotonic()-start


def certify(hist, h, etas):
    low, high = pi_intervals()
    last_bad = {str(eta):0 for eta in etas}
    minimum_j = max(h, 2)
    # At t>=1 each pi*t^2/4 - alpha*t used here is increasing.
    for start in range(minimum_j, len(hist)-1, 200000):
        stop = min(start+200000, len(hist)-1)
        j = np.arange(start, stop, dtype=np.int64)
        jl, jr = j, j+1
        left = low.numerator*jl*jl - 4*h*h*low.denominator*hist[jr]
        right = high.numerator*jr*jr - 4*h*h*high.denominator*hist[jl]
        for eta in etas:
            alpha, beta = 1-eta, 1+eta
            # Difference is formed BEFORE multiplying by the small eta denominator.
            # Verify that this final multiplication also fits int64.
            if max(int(np.max(np.abs(left))), int(np.max(np.abs(right))))*max(alpha.denominator,beta.denominator)>=2**63-1:
                raise OverflowError("corridor multiplication")
            if max(4*h*low.denominator*alpha.numerator*int(jl[-1]),
                   4*h*high.denominator*beta.numerator*int(jr[-1])) >= 2**63-1:
                raise OverflowError("corridor RHS multiplication")
            lower_fail = left*alpha.denominator < 4*h*low.denominator*alpha.numerator*jl
            upper_fail = right*beta.denominator > 4*h*high.denominator*beta.numerator*jr
            bad = np.flatnonzero(lower_fail | upper_fail)
            if len(bad):
                last_bad[str(eta)] = int(j[int(bad[-1])])+1
    return [{"eta":str(eta),
             "certified_start":(str(Q(max(minimum_j,last_bad[str(eta)]),h))
                                if last_bad[str(eta)]<len(hist)-1 else None),
             "certified_end":str(Q(len(hist)-1,h)),
             "largest_q_for_subtraction":str((1-eta)/(1+eta)),
             "tail_extension":"NOT included: use separately audited analytic discrepancy bound"}
            for eta in etas]


def paired_deficit_certificates(hist,h,qmax_values):
    """Compare outer lower d directly against the PREFIX maximum of inner d.

    All d-bounds are rounded OUTWARDS to integer multiples of 1/(4h^2).
    This retains local arithmetic cancellation lost by symmetric corridors.
    A successful outer bin certifies every 0<q<=qmax simultaneously.
    """
    low, high = pi_intervals()
    size = len(hist)-1
    upper = np.empty(size,dtype=np.int64)
    for start in range(0,size,200000):
        stop=min(start+200000,size)
        j=np.arange(start,stop,dtype=np.int64)
        num=high.numerator*(j+1)*(j+1)-4*h*h*high.denominator*hist[j]
        upper[start:stop]=-((-num)//high.denominator)
    np.maximum.accumulate(upper,out=upper)
    last_bad={str(q):0 for q in qmax_values}
    for start in range(h,size,200000):
        stop=min(start+200000,size)
        j=np.arange(start,stop,dtype=np.int64)
        num=low.numerator*j*j-4*h*h*low.denominator*hist[j+1]
        lower=num//low.denominator
        for q in qmax_values:
            if q.numerator*(int(j[-1])+1)>=2**63-1:
                raise OverflowError("prefix q product")
            inner_last=(q.numerator*(j+1))//q.denominator
            bad=np.flatnonzero(lower<upper[inner_last])
            if len(bad):
                last_bad[str(q)]=int(j[int(bad[-1])])+1
    return [{"qmax":str(q),
             "all_smaller_positive_q":True,
             "certified_start":(str(Q(max(h,last_bad[str(q)]),h)) if last_bad[str(q)]<size else None),
             "certified_end":str(Q(size,h)),
             "tail_extension":"NOT included: use separately audited analytic discrepancy bound"}
            for q in qmax_values]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--radius", type=int, default=2000)
    p.add_argument("--bins-per-unit", type=int, default=100)
    p.add_argument("--eta", action="append", default=[])
    p.add_argument("--output")
    p.add_argument("--paired",action="store_true")
    args = p.parse_args()
    etas = [Q(x) for x in (args.eta or ["1/10","9/100","2/25","7/100","3/50","1/20"])]
    if not all(Q(0)<x<Q(1,2) for x in etas):
        raise ValueError("eta must lie in (0,1/2)")
    hist, pairs, elapsed = build_histogram(args.radius,args.bins_per_unit)
    report={"status":"exact finite continuum corridors; analytic tails not included",
            "radius":args.radius,"bins_per_unit":args.bins_per_unit,
            "pairs":pairs,"seconds":elapsed,"histogram_bytes":hist.nbytes,
            "corridors":certify(hist,args.bins_per_unit,etas)}
    if args.paired:
        report["paired_deficit_certificates"]=paired_deficit_certificates(
            hist,args.bins_per_unit,[Q(4,5),Q(81,100),Q(41,50),Q(83,100),
                                    Q(21,25),Q(17,20),Q(43,50),Q(7,8),
                                    Q(9,10),Q(23,25),Q(19,20)])
    print(json.dumps(report,indent=2), flush=True)
    if args.output:
        Path(args.output).write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
