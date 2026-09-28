"""Independent last finite-gate replay, .90<=q<=.98, no generator imports.

Uses the previously independent replay's elementary assertion/floor
helpers, but constructs a fresh complete sorted lattice multiset and
rechecks every saved target index by a suffix minimum.
"""

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from time import monotonic
import json

import numpy as np

from replay_compact_band_independent import require, floor, PLO, PHI

CAP=2000000
SCALE=2**40


def main():
    started=monotonic()
    require(SCALE*CAP<2**63 and 64000**2*CAP<2**63,"integer array products")
    modes=np.asarray(sorted(m*m+n*n
                           for m in range(1,isqrt(CAP)+1)
                           for n in range(1,isqrt(CAP-m*m)+1)),dtype=np.int64)
    energies=np.arange(CAP+1,dtype=np.int64)
    counts=np.searchsorted(modes,energies,side="right")
    for e in (0,1,2,5,13,20,25,32,100,40401,250000,500000,1000000,CAP):
        direct=sum(isqrt(e-m*m) for m in range(1,isqrt(e)+1))
        require(int(counts[e])==direct,"direct row count at "+str(e))
    require(len(modes)==1569392,"complete ordered multiset")
    data=json.loads(Path(__file__).with_name("certificate_final.json").read_text(encoding="utf-8-sig"))
    records=data["records"]
    require(len(records)==1692,"interval count")
    previous=F(9,10)
    pairs=maximum=bypasses=0
    subtotals={"paired_tail":[0,0],"X_ge_8_tail":[0,0]}
    next_notice=monotonic()+20
    for number,record in enumerate(records,1):
        l,h=map(F,record["interval"])
        require(previous==l<h<=F(49,50),"gapless closed size cover")
        previous=h
        if h<=F(23,25):
            key="paired_tail"
            tail=F(163542,125)
            localcap=CAP
            kmax=floor(PHI*(1-l*l)*tail*tail/4)
            require(record["tail_gate"]=="paired-deficit tail","paired gate label")
        else:
            key="X_ge_8_tail"
            require(l>=F(23,25),"no interval straddles analytic gate")
            tail=8/(1-h)
            localcap=250000
            # For each actual q, the X=8 Weyl index threshold is
            # 16*pi*(1+q)/(1-q), increasing q. Use the RIGHT endpoint.
            kmax=floor(16*PHI*(1+h)/(1-h))
            require(record["tail_gate"]=="q>=23/25, X>=8 analytic cap tail (audit required)","analytic gate label")
        require(F(record["tail_upper"])==tail and record["cap"]==localcap,"tail expression and local cap")
        require(l.denominator<=64000,"reference denominator")
        mu=floor(SCALE*PLO*(1-h*h)/4)
        core=floor(PLO*(1+l)/(8*(1-l)))
        require(kmax==record["kmax"] and 2<=kmax<localcap,"finite cutoff")
        require(mu==record["mu"] and 0<mu<SCALE,"downward area coefficient")
        require(core==record["core_max"],"core bypass")
        tau=energies[1:localcap+1]
        inner=(l.numerator*l.numerator*tau)//(l.denominator*l.denominator)
        gap_index=counts[tau-1]-counts[inner]+1
        require(bool(np.all(gap_index>=1)),"positive gap indices")
        suffix_best=np.minimum.accumulate(gap_index[::-1])[::-1]
        ks=np.arange(1,kmax+1,dtype=np.int64)
        thresholds=(SCALE*ks+mu-1)//mu
        require(bool(np.all((thresholds>=1)&(thresholds<=localcap))),"candidate thresholds inside census")
        passed=(ks<=max(2,core))|(suffix_best[thresholds-1]<=ks)
        require(ks[~passed].tolist()==record["remaining"]==[],"every target index passes")
        pairs+=kmax
        maximum=max(maximum,kmax)
        bypasses+=max(2,core)
        subtotals[key][0]+=1
        subtotals[key][1]+=kmax
        if monotonic()>=next_notice:
            print(json.dumps({"phase":"independent final suffix replay","intervals":number,"pairs":pairs,"seconds":round(monotonic()-started,3)}),flush=True)
            next_notice=monotonic()+20
    require(previous==F(49,50),"final endpoint")
    require(pairs==241379723 and data["pairs"]==pairs,"complete finite pair count")
    require(data["intervals"]==len(records),"saved interval total")
    require(data["remaining_indices"]==[] and data["remaining"]==[] and data["remaining_interval_count"]==0,"saved empty residual summaries")
    print(json.dumps({"status":"PASS_INDEPENDENT_FINAL_FINITE_REPLAY",
                      "intervals":len(records),"pairs":pairs,"maximum_index":maximum,
                      "ordered_lattice_modes":len(modes),"core_or_universal_bypass_pairs":bypasses,
                      "gate_interval_pair_totals":subtotals,"remaining":[],
                      "seconds":round(monotonic()-started,3),
                      "analytic_gate_status":"X>=8 theorem requires separate audit; finite replay does not assert it"},indent=2))


if __name__=="__main__":
    main()
