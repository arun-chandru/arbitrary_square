"""Independent small-case audit; does not repeat the R=50000 census."""

from fractions import Fraction as Q
from math import ceil, floor, isqrt

from certify_square_corridor_histogram import build_histogram, certify, pi_intervals, paired_deficit_certificates


def ordered_histogram(R,h):
    bins=[0]*(R*h+1)
    for m in range(1,R+1):
        for n in range(1,R+1):
            if m*m+n*n<=R*R:
                square=(m*m+n*n)*h*h
                index=isqrt(square)
                if index*index<square:
                    index+=1
                bins[index]+=1
    count=0
    for j in range(len(bins)):
        count+=bins[j]
        bins[j]=count
    return bins


def independent_starts(hist,h,etas):
    pi_low,pi_high=pi_intervals()
    starts={str(eta):Q(max(h,2),h) for eta in etas}
    last_bad={str(eta):0 for eta in etas}
    for j in range(max(h,2),len(hist)-1):
        left,right=Q(j,h),Q(j+1,h)
        for eta in etas:
            lower=pi_low*left*left/4-hist[j+1]-(1-eta)*left
            upper=pi_high*right*right/4-hist[j]-(1+eta)*right
            if lower<0 or upper>0:
                starts[str(eta)]=right
                last_bad[str(eta)]=j+1
    endpoint=Q(len(hist)-1,h)
    return {eta:None if last_bad[eta]==len(hist)-1 else value for eta,value in starts.items()}


def independent_paired_starts(hist,h,qvalues):
    pi_low,pi_high=pi_intervals()
    upper=[]
    for j in range(len(hist)-1):
        u=ceil(4*h*h*(pi_high*Q(j+1,h)**2/4-hist[j]))
        upper.append(max(u,upper[-1] if upper else u))
    starts={str(q):Q(1) for q in qvalues}
    for j in range(h,len(hist)-1):
        lower=floor(4*h*h*(pi_low*Q(j,h)**2/4-hist[j+1]))
        for q in qvalues:
            if lower<upper[floor(q*(j+1))]:
                starts[str(q)]=Q(j+1,h)
    endpoint=Q(len(hist)-1,h)
    return {q:None if value==endpoint else value for q,value in starts.items()}


def main():
    cases=0
    bins=0
    etas=[Q(1,10),Q(9,100),Q(2,25),Q(1,20)]
    for R in (2,3,17,61):
        for h in (1,3,7,100):
            actual,_,_=build_histogram(R,h)
            expected=ordered_histogram(R,h)
            if list(actual)!=expected:
                raise AssertionError((R,h,"histogram"))
            starts=independent_starts(expected,h,etas)
            for row in certify(actual,h,etas):
                actual_start=None if row["certified_start"] is None else Q(row["certified_start"])
                if actual_start!=starts[row["eta"]]:
                    raise AssertionError((R,h,row,starts))
            qvalues=[Q(4,5),Q(41,50),Q(19,20)]
            paired_starts=independent_paired_starts(expected,h,qvalues)
            for row in paired_deficit_certificates(actual,h,qvalues):
                actual_start=None if row["certified_start"] is None else Q(row["certified_start"])
                if actual_start!=paired_starts[row["qmax"]]:
                    raise AssertionError((R,h,row,paired_starts))
            cases+=1
            bins+=len(expected)
    print(f"PASS: {cases} complete small histograms, {bins} bins; exact Fraction corridor AND paired-prefix inequalities agree.")
    # Separately exhibit a generic reporting issue, not affecting the saved
    # 50000 run: if the final bin fails, a degenerate endpoint is not certified.
    hist,_,_=build_histogram(2,10)
    row=certify(hist,10,[Q(1,100)])[0]
    if row["certified_start"] is not None or row["certified_end"]!="2":
        raise AssertionError("degenerate-report fix failed")
    # N(2)=1, so d(2)=pi-1>2.02 since pi>3.02.
    if not pi_intervals()[0]>Q(302,100):
        raise AssertionError("pi check")
    print("DEGENERATE REPORT FIX CONFIRMED: R=2,h=10,eta=.01 now correctly reports no certified interval.")


if __name__=="__main__":
    main()
