"""Independent exact replay for .8<=q<=.83; imports no certificate generator.

Sorted ordered lattice modes + suffix minima replace the generator's
cumulative-max frontier. Every saved index is checked. Only the separately
repaired k=5 band may remain. No files are written; safe under -B -O.
"""

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

import numpy as np


CAP=250000
SCALE=2**40
PLO,PHI=F(103993,33102),F(104348,33215)
TAIL=F(977,4)
REPAIR=(F(643,800),F(129,160))


def require(test,message):
    if not test:
        raise RuntimeError(message)


def floor(q):
    return q.numerator//q.denominator


def certify_pi():
    def atan(x):
        s=sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(12)),F(0))
        return s,s+x**25/25
    a,b=atan(F(1,5)),atan(F(1,239))
    require(PLO<16*a[0]-4*b[1]<16*a[1]-4*b[0]<PHI,"independent Machin pi bracket")


def main():
    certify_pi()
    require(SCALE*CAP<2**63 and 64000**2*CAP<2**63,"integer array range")
    modes=np.asarray(sorted(m*m+n*n
                           for m in range(1,isqrt(CAP)+1)
                           for n in range(1,isqrt(CAP-m*m)+1)),dtype=np.int64)
    require(len(modes)==195837,"complete ordered disk census")
    energies=np.arange(CAP+1,dtype=np.int64)
    counts=np.searchsorted(modes,energies,side="right")
    for e in (0,1,2,5,18,25,32,40,40401,100000,CAP):
        direct=sum(isqrt(e-m*m) for m in range(1,isqrt(e)+1))
        require(int(counts[e])==direct,"independent direct row count at "+str(e))
    tau=energies[1:]
    data=json.loads(Path(__file__).with_name("certificate_upper_midband.json").read_text(encoding="utf-8-sig"))
    records=data["records"]
    require(F(data["tail"])==TAIL and len(records)==120,"tail and interval count")
    previous=F(4,5)
    pairs=maximum=core_bypasses=0
    unmatched=[]
    repair_spans=[]
    for record in records:
        l,h=map(F,record["interval"])
        require(previous==l<h<=F(83,100),"closed interval coverage")
        previous=h
        require(l.denominator<=64000 and l.numerator<l.denominator,"reference arithmetic range")
        kmax=floor(PHI*(1-l*l)*TAIL**2/4)
        mu=floor(SCALE*PLO*(1-h*h)/4)
        # A(q)*lambda_1/(4*pi) >= pi*(1+q)/(8*(1-q)), increasing q.
        core_max=floor(PLO*(1+l)/(8*(1-l)))
        require(kmax==record["kmax"] and 2<=kmax<CAP,"finite index cutoff")
        require(mu==record["mu"] and 0<mu<SCALE,"downward area coefficient")
        require(core_max==record["core_max"],"common-core bypass")
        inner=(l.numerator*l.numerator*tau)//(l.denominator*l.denominator)
        gap_index=counts[tau-1]-counts[inner]+1
        require(bool(np.all(gap_index>=1)),"positive gap indices")
        suffix_best=np.minimum.accumulate(gap_index[::-1])[::-1]
        ks=np.arange(1,kmax+1,dtype=np.int64)
        required_tau=(SCALE*ks+mu-1)//mu
        require(bool(np.all((required_tau>=1)&(required_tau<=CAP))),"thresholds inside census")
        passed=(ks<=max(2,core_max))|(suffix_best[required_tau-1]<=ks)
        missing=ks[~passed].tolist()
        require(missing==record["remaining"],"independent missing-index list")
        if missing:
            require(missing==[5],"only k=5 needs repair")
            require(REPAIR[0]<=l<h<=REPAIR[1],"entire residual interval inside repair")
            repair_spans.append((l,h))
            unmatched.append(record)
        pairs+=kmax
        maximum=max(maximum,kmax)
        core_bypasses+=max(2,core_max)
    require(previous==F(83,100),"final endpoint")
    require(pairs==1888608,"every index/interval pair")
    require(data["intervals"]==len(records) and data["pairs"]==pairs,"saved totals")
    require(data["remaining_indices"]==[5] and data["remaining"]==unmatched,"saved residual summary")
    require(repair_spans[0][0]==REPAIR[0] and repair_spans[-1][1]==REPAIR[1],"repair endpoints")
    for first,second in zip(repair_spans,repair_spans[1:]):
        require(first[1]==second[0],"contiguous residual band")
    print(json.dumps({"status":"PASS_INDEPENDENT_UPPER_MIDBAND_REPLAY",
                      "intervals":len(records),"pairs":pairs,"maximum_index":maximum,
                      "core_or_universal_bypass_pairs":core_bypasses,
                      "residual_index":5,"residual_interval_pairs":len(repair_spans),
                      "exact_residual_band":[str(x) for x in REPAIR],
                      "method":"sorted ordered lattice modes; inclusive search; suffix-minimum gap-index test",
                      "dependency":"separate all-angle k5 repair and paired-deficit tail at t>=977/4"},indent=2))


if __name__=="__main__":
    main()
