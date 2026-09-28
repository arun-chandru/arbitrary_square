# Pólya's inequality for a square with an arbitrary square hole

R. Arun Chandru · 28 September 2026

## Extended abstract

Imagine a square membrane held fixed at its outer boundary and at the boundary
of a square hole cut out of its interior. Its vibration frequencies are
determined by the Dirichlet Laplacian. Weyl's law predicts their average density
at high frequencies. Pólya's conjecture makes a stronger assertion: the number
of frequencies below any specified energy should never exceed the leading
area-based prediction.

This paper establishes that inequality for a square with **any square hole
contained in its closure**, provided the remaining area is positive. The hole
need not be centered or parallel to the outer square. The result applies to
every eigenvalue, counting multiplicity,
not just to sufficiently high frequencies. No positive uniform lower bound
on the hole's distance from the outer boundary is imposed. Boundary contact
and disconnected remaining domains are included by a final one-sided
approximation using domain monotonicity and convergence of the areas.

The proof combines three kinds of information. First, the explicitly known
spectra of the outer square and isolated square hole give comparisons that
are independent of the hole's position. Second, information from a few low
square modes and carefully controlled changes of the actual domain treats
the finite range where those comparisons alone are insufficient. Third,
when the hole nearly fills the square, a decomposition into annular sectors
and corners retains the area and boundary information needed for a uniform
estimate. Two conformal corner models are used with their changed mass
weights and endpoint conditions explicitly accounted for.

The unbounded parameter ranges are treated analytically. The remaining
bounded obligations are exact integer and rational inequalities, verified
on whole parameter intervals. These computations do not approximate
eigenvalues of the perforated square and do not infer continuum statements
from sampled shapes. The accompanying programs expose the finite tests,
their inputs, and the arithmetic checks needed to reproduce them.

Two further results are included. A sharp threshold is obtained for a
particular complementary-count criterion for multiple square holes: the
sum of their side lengths may be at most
`sqrt(1 - 12/(5*pi))` times the outer side. The threshold is sharp for this
scalar criterion, not asserted to be an optimal geometric restriction.
These square holes may touch one another or the outer boundary, provided
their interiors are pairwise disjoint. A separate argument treats an
arbitrarily displaced, strictly interior parallel rectangular hole when
both side ratios are at least `11/12`.




## Compiling the paper

`r_arun_chandru_Pólyas_inequality_square_with_arbitrary_square_hole.pdf` is the complete paper, including its appendices and references.
`manuscript.tex` is its **single, self-contained LaTeX source**. 

With a standard TeX Live or MiKTeX installation, run:

```text
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

Alternatively, run `tectonic manuscript.tex`. The source uses `amsart`, AMS math
packages, Latin Modern, `geometry`, `booktabs`, `array`, `longtable`, TikZ
(`calc` and `arrows.meta`), `microtype`, `xurl`, and `hyperref`. No BibTeX or
external figure-generation command is needed. The Python computations are
separate from compilation: compiling the paper does not rerun them.


## Contents of the verification supplement

The executable source is deliberately kept separate from the LaTeX source.
The five files below contain rational witness inputs; a successful status
field in an input is never accepted in place of recomputing the tests.

| Exact witness input | Independent replay |
|---|---|
| `all_angles_certificate.json` | `replay_all_angles_band.py` |
| `certificate_midband.json` | `replay_midband_independent.py` |
| `certificate_upper_midband.json` | `replay_upper_midband_independent.py` |
| `certificate_compact.json` | `replay_compact_band_independent.py` |
| `certificate_final.json` | `replay_final_compact_independent.py` |

These reconstruct square-lattice counts, check the complete closed size
cover, recompute finite-index cutoffs, and check all interval/index
obligations. The ninth-, eleventh-, and fifth-mode residual bands are
handled by the explicitly proved angular comparisons, not erased from the
input records.

The other source files have the following roles:

| Programs | Role |
|---|---|
| `verify_deficits.py`, `verify_corridor.py` | Multiple-hole scalar comparisons and the radius-201 deficit corridor. |
| `verify_rectangular.py` | Five finite scalar comparisons for the rectangular-hole conclusion and the midpoint estimate. |
| `verify_rotation_gram.py` | Exact low-mode Gram bound on 176 closed intervals. |
| `verify_low_mode_rotation_repairs.py`, `verify_fifth_mode_upper_repair.py` | Reconstruct the exact angular covers. Their generated leaf lists are outputs, not required inputs. |
| `verify_radial_scalar.py` | The scalar radial deficit on 75 closed intervals, with its analytic-tail constant. |
| `verify_cap_closure_constants.py`, `verify_cap_q92_x8_constants.py` | Rational constants and closing margins for the two cap estimates. These do not replace the written operator proofs. |
| `verify_circle_tail_constants.py` | Rational constants for the explicit circle estimate and analytic continuations. |
| `certify_square_corridor_histogram.py` | Reconstruct the three large, exact paired-deficit histograms. |
| `audit_square_corridor_histogram_small.py` | Independently check the histogram construction and rounding conventions at small radii. |
| `verify_all.py` | Sequential execution, failure reporting, and reproducibility logs. |

The three supplied `square_corridor_paired_*.json` reports are compact
reference outputs from completed computations, not lattice-count inputs.
Reproduction rebuilds the counts instead of trusting the reports; the runner
then requires the fresh mathematical report fields to match these reference
outputs, excluding the timing field. Files in `results/` record the
release replay; they are outputs. Timing and runtime metadata need not be
identical on another computer. `SHA256SUMS.txt` identifies the distributed
files; it is an integrity record, not a mathematical certificate.

## Reproducing the checks

Use Python 3.10 or newer and NumPy. All other imports are from the standard
library. No network access, symbolic algebra package, numerical eigensolver,
Bessel-function implementation, or external research folder is required.
All source files should remain together, as a few import elementary helpers
from a sibling verifier. From this directory run:

```text
python -B -O verify_all.py --group all
```

The runner executes sequentially and writes logs and reports to `results/`
by default. To preserve the supplied release reports, select another output
directory:

```text
python -B -O verify_all.py --group all --output-root my_results
```

Individual groups are `arithmetic`, `geometry`, `inherited`, and `histograms`.
The `--group` option can be repeated. Omitting the large histogram runs is
useful for a partial check, but is not a complete reproduction.
Each verifier can also be run directly. Acceptance checks raise explicit
exceptions and remain active under Python's `-O` flag.

The three histogram commands, if run directly, are:

```text
python -B -O certify_square_corridor_histogram.py --radius 50000 --bins-per-unit 100 --paired --eta 1/14 --eta 1/19 --output replay_50000.json
python -B -O certify_square_corridor_histogram.py --radius 5000 --bins-per-unit 1000 --paired --eta 1/19 --eta 1/24 --output replay_5000.json
python -B -O certify_square_corridor_histogram.py --radius 90000 --bins-per-unit 100 --paired --eta 1/24 --output replay_90000.json
```

Integer and rational acceptance bounds, including signed-64-bit safety
checks, are part of the algorithms. Where floating-point seeds accelerate
integer bin construction, the integer inequalities defining the resulting
bins are checked explicitly; an unchecked floating approximation is not
an acceptance criterion.

## Computational size and limits of verification

The five finite-band replays cover 3,635 closed size intervals and
529,530,599 interval/index obligations, including the stated analytic
low-mode alternatives. Their largest local squared-radius cap is 2,000,000.
The three histograms process respectively 981,740,308, 9,816,727, and
3,180,849,242 unordered positive lattice pairs. This is finite exhaustive
arithmetic, not a search over approximate domain spectra.

The largest histogram array itself uses 72,000,008 bytes; prefix arrays and
temporary arrays require additional memory. The implementation uses
vectorized rows and arrays of order the radius times the bin density, not
a full radius-squared array. Run the large jobs sequentially. Runtime depends
on the machine and NumPy build; the release logs give the actual run times
on the verification machine, not a hardware-independent guarantee.

The release was checked with Python 3.12.14 and NumPy 2.3.5. All 19
sequential checks passed in 260.156 seconds on that machine. The three
full histogram replays reproduced the supplied mathematical report fields
exactly, excluding elapsed time. The complete execution record is
`results/verification_results.json`.

