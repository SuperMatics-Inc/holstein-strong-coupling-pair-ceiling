# A strong-coupling ceiling of 6000 t^2/U for electronic pair correlations in the two-dimensional Holstein model

Davis Rens, SuperMatics Inc. ([ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)).
Version 1.0.0.

Paper DOI: <https://doi.org/10.5281/zenodo.24414616>
Replay supplement DOI: <https://doi.org/10.5281/zenodo.24414618>
Source: <https://github.com/SuperMatics-Inc/holstein-strong-coupling-pair-ceiling>

## What this deposit is

This repository carries the preprint `paper/main.pdf` (19 pages) and the complete replay supplement for it.
The paper proves an exponential clustering bound for the equal-time on-site electronic pair correlations of
the spinful Holstein model on the two-dimensional square torus, with a single hypothesis on `beta t^2/U`, no
restriction on `t/U`, and a bound uniform in the phonon frequency, in the torus side and in the chemical
potential.

## What is proved

Write `U = F^2/K > 0` and let `beta` be the inverse temperature. There is no hypothesis on the ratio `t/U`, and
every statement is uniform in `Omega > 0`, in the side `N >= 4` and in `mu` real.

1. For every `N >= 4`, if `beta t^2/U <= 1/6000` then every pair `x != y` of the torus satisfies
   `|<Delta_x^dagger Delta_y>| <= (57/500)(99/100)^d` in the torus distance `d >= 1`. The rows of the
   correlation matrix are then uniformly summable, bounded by `1 + 4Aq/(1-q)^2 <= 4515.4` with
   `A = 57/500` and `q = 99/100`, and the clustering temperature of the definition in the paper obeys
   `k_B T_pair <= 6000 t^2/U`.
2. If `N` is even the same decay holds under the weaker hypothesis `beta t^2/U <= 1/4500`, hence
   `k_B T_pair <= 4500 t^2/U` on even tori.

The statement also holds in canonical thermodynamic limits at interior density. It imposes no smallness
condition on `t/U` and does not assume the adiabatic limit: the bound is unchanged as `Omega -> 0`.

## Where the constants come from

The mechanism is exact at every step, and each step is a separate section of the paper:

1. an exact signless majorant of the oscillator contraction, in the Lang-Firsov form, positive in the dipole
   kernel and therefore valid for every phonon frequency;
2. a connected-support expansion of the pair correlation with a signless trajectory majorant;
3. the exact sixteen-state single-bond spectrum, with its sector multiplicities and the half-filling
   supremum of the bond trace;
4. a bridge parity that doubles the effective decay factor, an explicit four- or five-edge colouring of the
   connected supports, and a Schatten Holder decoupling that costs a factor `c <= 5`;
5. exact enumeration of connected supports with at most six edges and a four-ary Raney tree tail;
6. the Kotecky-Preiss vertex criterion at `B = 557/500`.

## What this adds, and what it does not claim

Against the previously released result for the same model (`10.5281/zenodo.24380128`), which gave
`k_B T_pair <= 2^27 t^2/U = 134 217 728 t^2/U` with prefactor `240/59` and ratio `5/64`: the ceiling improves
by a factor `22 369` on all tori and `29 826` on even tori, the prefactor by a factor `35.68`, and the
exponential rate `|ln q|` by a factor `254`, at the price of a much slower decay ratio, `99/100` against
`5/64`.

Against the dimension-tracking preprint (`10.5281/zenodo.24412193`), whose two-dimensional constant is
`15 925 248` with ratio `5/24` and prefactor `360/19`: the ceiling here is smaller by `2654` on all tori and
`3539` on even tori, the prefactor by `166`, and the exponential rate by `156`. That preprint treats the
general cubic case by a different route and gives `63 700 992` in three dimensions; a two-dimensional
argument of the present kind does not transfer without redoing the colouring and the tree count, and no claim
in three dimensions is made here.

The limitations are stated in the paper and hold here as well.

- The bound is on a clustering temperature of one correlator. It is not the critical temperature of
  anything, no phase transition is proved or excluded, and nothing shows that the model is superconducting at
  any temperature.
- The decay rate is weak, `|ln q| = 0.01005`, so the bound is a statement about uniform summability of the
  correlation matrix and about an exact threshold in `beta t^2/U`, not a quantitative decay estimate at
  moderate distance.
- The proof covers the two-dimensional square lattice. The colouring, the small-edge enumeration and the
  four-ary tree count are specific to the square grid.
- The constants `6000` and `4500` are not optimal. Optimising the continuous parameters of the same scheme
  with the combinatorial tables held fixed puts its floor at about `5314` on all tori and `4000` on even
  tori, and the certified failures at `C = 5800` and `C = 4300` are failures of the tail test, evidence about
  the scheme and not about the model.
- The certificate program proves arithmetic, not analysis.

## Contents

```
MANIFEST.sha256          SHA-256 of every tracked file except this one
CITATION.cff             citation metadata for the preprint and this supplement
LICENSE                  MIT License, for the code
paper/main.pdf           the reviewed manuscript, released unmodified
paper/main.tex           LaTeX source of that manuscript
paper/refs.bib           bibliography
replay/replay.py         the standalone certificate program
replay/README.md         what the program verifies and how to run it
replay/MANIFEST.sha256   SHA-256 of replay.py and replay/README.md, pinned since review
```

`paper/main.pdf` is the reviewed artifact and is released byte for byte; its SHA-256 digest is
`bf0b5b56805fa8d2f0c2f806d962fa3b26fc5ea0668445591b28a78973897fab` and `paper/main.tex` is
`cf7576c78a115caec8e6280b096d9a0ca8669ef2373eee7020fc41b2332021a2`. Continuous integration verifies those
bytes, the text they carry and the replay program; it does not rebuild the PDF from `main.tex`, so the deposit
is exactly the reviewed document.

## Replay

`replay/replay.py` is one standalone program written against the Python standard library only: no data file,
no network access, no third-party package and no floating-point proof step. Python 3.13 is what the review
used; any Python 3.9 or newer should run it.

```
python replay/replay.py
python replay/replay.py --mutants
```

The default run prints seven certified lines and exits zero only if all of them hold, from `16STATE_HALF_FILLING`
through `WEIGHTED_2AB_CURVE`: the sixteen-state CAR and bond-trace identities with the two half-filling
suprema; the exact torus enumeration at `k <= 6`, the four-ary coefficients `R4` and the two-variable Raney
tree check; the 102 certified outward Fraction intervals of the pair bounds over 22 class slots; the budgets
of the ceiling at `C = 6000` on all tori and `C = 4500` on even tori, with the canonical source transfer; the
ten points of the unconditional trade-off curve; and the five certified points of the weighted two-antibond
route. The `--mutants` run repeats the checks in fresh child processes with one guarded claim corrupted at a
time, twelve times, and ends with `CERTIFIED: all 12 mutants rejected with exit 1`; two of the corrupted
claims are the ceilings `C = 5800` on all tori and `C = 4300` on even tori, which the same scheme must refuse.

```
$ python replay/replay.py
16STATE_HALF_FILLING: S_JW*T_e*S_JW==B_e, 5-sector traces & 2 supremum rational identities verified
ENUM_AND_RANEY: T_4..T_7 exact k<=6 counts dominated; R4=[1,4,18,88,455,2448,13566]; 2-var Raney (kp+ks<=4) verified
PE_INTERVALS: 22 class bounds certified across 102 adaptive outward Fraction intervals
MAIN_ALL_TORI: eps0=inf, C=6000, q_pair=99/100, A_pair=57/500, S_phys<=0.104410, S_src<=0.000857, S_tot<=0.105268 <= lnB>=0.107957 (22369.6x tighter than the earlier 2D release C=2^27; 45812984x tighter than the earlier conditional release C=2^38)
MAIN_EVEN_TORI: eps0=inf, C=4500, q_pair=99/100, A_pair=57/500, S_phys<=0.101996, S_src<=0.000838, S_tot<=0.102834 <= lnB>=0.107957 (29826.2x tighter than the earlier 2D release C=2^27; 61083979x tighter than the earlier conditional release C=2^38)
UNCONDITIONAL_CURVE: 10 rational (mode, q_pair, C) points certified
WEIGHTED_2AB_CURVE: 5 rational points certified (C=85000 at eps0=173/1000000 up to C=765000 at eps0=1537/1000000)
```

`replay/README.md` is the guide pinned since review: it records the commands, the observed exit codes and
what each line certifies. The program proves arithmetic, not analysis.

`paper/main.pdf` was built from `paper/main.tex` with the toolchain recorded in the review: TeXLive (2026),
Python 3.13.9, `uv` 0.9.5, `sympy` 1.14.0 and `mpmath` 1.3.0. The deposit releases the built PDF; it does not
rebuild it.

## Licences

The code in this repository, that is `replay/replay.py` and the workflow, is released under the MIT License
(`LICENSE`). The manuscript text, `paper/main.pdf`, `paper/main.tex` and `paper/refs.bib`, is released under
the Creative Commons Attribution 4.0 International license (CC-BY-4.0), as on the paper record. The replay
supplement record is released under the MIT License.

## Status and disclosures

AI agents generated the proofs, code and drafts under D. Rens' direction, and separate AI systems carried out
the independent reviews of the manuscript and of the certificate program. This preprint has not been peer
reviewed.
