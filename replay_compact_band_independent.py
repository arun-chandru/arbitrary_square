"""Independent all-pair replay for .83<=q<=.9, using sorted lattice modes.

Imports neither generator nor frontier code. The certificate's future-
energy maximum is checked instead as a suffix minimum of gap indices.
One process, bounded O(cap) arrays, no files written, safe under -B -O.
"""

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from time import monotonic
import json

import numpy as np

CAP=2000000
SCALE=2**40
PLO,PHI=F(103993,33102),F(104348,33215)
BANDS=[(F(83,100),F(21,25),F(977,4),250000),
       (F(21,25),F(43,50),F(38427,100),250000),
       (F(43,50),F(7,8),F(2317,4),500000),
       (F(7,8),F(9,10),F(67111,50),2000000)]


def require(test,message):
    if not test:
        raise RuntimeError(message)


def floor(x):
    return x.numerator//x.denominator


def main():
    started=monotonic()
    def atan(x):
        s=sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(12)),F(0))
        return s,s+x**25/25
    a,b=atan(F(1,5)),atan(F(1,239))
    require(PLO<16*a[0]-4*b[1]<16*a[1]-4*b[0]<PHI,"independent pi bracket")
    require(SCALE*CAP<2**63 and 64000**2*CAP<2**63,"integer product ranges")
    modes=np.asarray(sorted(m*m+n*n
                           for m in range(1,isqrt(CAP)+1)
                           for n in range(1,isqrt(CAP-m*m)+1)),dtype=np.int64)
    energies=np.arange(CAP+1,dtype=np.int64)
    counts=np.searchsorted(modes,energies,side="right")
    for e in (0,1,2,5,13,20,25,32,100,40401,250000,500000,1000000,CAP):
        direct=sum(isqrt(e-m*m) for m in range(1,isqrt(e)+1))
        require(int(counts[e])==direct,"independent direct row count at "+str(e))
    require(len(modes)==int(counts[CAP]),"complete ordered multiset")
    data=json.loads(Path(__file__).with_name("certificate_compact.json").read_text(encoding="utf-8-sig"))
    records=data["records"]
    require(len(records)==1241,"interval total")
    previous=F(83,100)
    pairs=maximum=bypasses=0
    next_notice=monotonic()+20
    band_counts={str(a)+".."+str(b):[0,0] for a,b,_,_ in BANDS}
    for number,record in enumerate(records,1):
        l,h=map(F,record["interval"])
        require(previous==l<h<=F(9,10),"gapless closed interval coverage")
        previous=h
        matched=[x for x in BANDS if x[0]<=l<h<=x[1]]
        require(len(matched)==1,"one prescribed tail band")
        blo,bhi,tail,localcap=matched[0]
        require(F(record["tail"])==tail and record["cap"]==localcap,"prescribed tail and census cap")
        require(l.denominator<=64000,"reference arithmetic bound")
        kmax=floor(PHI*(1-l*l)*tail**2/4)
        mu=floor(SCALE*PLO*(1-h*h)/4)
        core=floor(PLO*(1+l)/(8*(1-l)))
        require(kmax==record["kmax"] and 2<=kmax<localcap,"finite index cutoff")
        require(mu==record["mu"] and 0<mu<SCALE,"downward area coefficient")
        require(core==record["core_max"],"common-core bypass")
        tau=energies[1:localcap+1]
        inner=(l.numerator*l.numerator*tau)//(l.denominator*l.denominator)
        gap_index=counts[tau-1]-counts[inner]+1
        require(bool(np.all(gap_index>=1)),"positive gap indices")
        suffix_best=np.minimum.accumulate(gap_index[::-1])[::-1]
        ks=np.arange(1,kmax+1,dtype=np.int64)
        thresholds=(SCALE*ks+mu-1)//mu
        require(bool(np.all((thresholds>=1)&(thresholds<=localcap))),"all thresholds inside local census")
        passes=(ks<=max(2,core))|(suffix_best[thresholds-1]<=ks)
        missing=ks[~passes].tolist()
        require(missing==record["remaining"]==[],"no missing index")
        pairs+=kmax
        maximum=max(maximum,kmax)
        bypasses+=max(2,core)
        bkey=str(blo)+".."+str(bhi)
        band_counts[bkey][0]+=1
        band_counts[bkey][1]+=kmax
        if monotonic()>=next_notice:
            print(json.dumps({"phase":"independent suffix replay","intervals":number,"pairs":pairs,"seconds":round(monotonic()-started,3)}),flush=True)
            next_notice=monotonic()+20
    require(previous==F(9,10),"final endpoint")
    require(pairs==261974004 and data["pairs"]==pairs,"all finite pairs replayed")
    require(data["intervals"]==len(records),"saved interval total")
    require(data["remaining_indices"]==[] and data["remaining"]==[] and data["remaining_interval_count"]==0,"saved residual summary")
    print(json.dumps({"status":"PASS_INDEPENDENT_COMPACT_BAND_REPLAY",
                      "intervals":len(records),"pairs":pairs,"maximum_index":maximum,
                      "ordered_lattice_modes":len(modes),"core_or_universal_bypass_pairs":bypasses,
                      "bands_interval_pair_totals":band_counts,
                      "remaining":[],"seconds":round(monotonic()-started,3),
                      "method":"sorted ordered lattice multiset, inclusive search, independent suffix minimum",
                      "tail_dependency":"Paired-deficit certificates on whole radius intervals, manuscript.tex"},indent=2))


if __name__=="__main__":
    main()
