# Supplement: exact certificate for the strong-coupling pair-clustering ceiling

This directory contains one certificate program, `replay.py`, that recomputes, in exact
rational and outward-interval arithmetic, every finite witness quoted in the preprint
"A strong-coupling ceiling of $6000\,t^{2}/U$ for electronic pair correlations in the
two-dimensional Holstein model" (Davis Rens, SuperMatics Inc., 2026).

The program certifies finite statements only: exact rational enclosures, exact operator
identities on small Fock spaces, exact lattice-animal counts, and exact checks of the
Kotecky-Preiss budgets. It does not certify the analytical many-body theorems; those are
proved in the preprint, and the program checks the arithmetic constants that appear in
their statements.

## Requirements

* `uv` 0.9.5 (any recent `uv` works)
* Python 3.13.9 (fetched and used by `uv`)
* `sympy` 1.14.0
* `mpmath` 1.3.0

Both packages are public on the Python Package Index and are installed by `uv` into a
throwaway environment. The program itself uses only the Python standard library, no
floating-point proof step, no network access and no stored data.

## Commands

Default run: every exact witness must certify.

```
uv run --no-project --with sympy --with mpmath python replay.py
```

Guarded mutant run: twelve deliberately corrupted variants of the certificate must each
be rejected, and the program must exit zero only if all twelve are rejected.

```
uv run --no-project --with sympy --with mpmath python replay.py --mutants
```

Observed exit codes, in a clean copy of this directory on Windows 11 with `uv` 0.9.5:

| command | exit code | output summary |
| --- | --- | --- |
| default | 0 | 7 witness lines, `16STATE_HALF_FILLING` to `WEIGHTED_2AB_CURVE` |
| `--mutants` | 0 | 12 `REJECTED <mutant>: <guard>` lines, then `CERTIFIED: all 12 mutants rejected with exit 1` |

A nonzero exit code means a certificate failed; in the mutant run it also means a
corrupted certificate was accepted, which is the failure mode that run is designed to
detect. Each mutant is run in a fresh child process, so an accepted mutant cannot
contaminate the evidence of the others.

## What is certified

Default run, in order:

1. `16STATE_HALF_FILLING`. On the sixteen-dimensional four-orbital Fock space of a single
   lattice edge, the Jordan-Wigner gauge transforms the physical hopping into the
   signless hopping as the exact operator identity `S_JW T_e S_JW = B_e`, together with
   the five sector traces of the bond weight and the rational identities used for the
   half-filling supremum of the bond transfer ratio.
2. `ENUM_AND_RANEY`. The exact counts of connected hopping supports with at most six
   edges on the tori `T_4` to `T_7`, their domination by the corresponding values of the
   four-ary tree tail `R4 = [1, 4, 18, 88, 455, 2448, 13566]`, and the two-variable
   (k-claw and k-star) Raney multigraph counts for `kp + ks <= 4` that bound the tail.
3. `PE_INTERVALS`. Twenty-two class bounds on the bond-decoupled ratio `P_e(c,y)` over
   102 adaptive outward `Fraction` intervals, each interval verified with rational
   Taylor and square-root enclosures rather than floating point.
4. `MAIN_ALL_TORI`. The Kotecky-Preiss budget at `eps0 = infinity`, `C = 6000`,
   `q = 99/100`, `A = 57/500` on all tori `N >= 4`: the physical penalty is at most
   `0.104410`, the canonical-source penalty at most `0.000857`, their sum at most
   `0.105268`, below the available radius `ln B >= 0.107957`. The line also prints the
   exact improvement factors against the earlier constant `C = 2^27` and the earlier
   conditional constant `C = 2^38`.
5. `MAIN_EVEN_TORI`. The same budget at `C = 4500` on tori of even side, with the sum of
   penalties at most `0.102834`.
6. `UNCONDITIONAL_CURVE`. The ten rational points of the unconditional trade-off curve
   between the ceiling `C` and the decay ratio `q_pair`, each certified with its own
   budget.
7. `WEIGHTED_2AB_CURVE`. The five rational points of the weighted two-`ab` obstruction,
   from `C = 85000` at `eps0 = 173/1000000` to `C = 765000` at `eps0 = 1537/1000000`.

The mutant run corrupts one guarded claim at a time: the bridge multiplicity `m = k`,
the colour count on three-edge claws, the colour count on odd tori, one interval in the
`P_e` block, the sign of the half-filling term, the parity classes on odd tori, the
ceiling lowered to `C = 5800` and to `C = 4300`, the canonical-source scale, the pair
prefactor, an endpoint of the weighted obstruction, and the ternary coefficient of the
Raney tail. Each corruption must be rejected.

## Manifest

`MANIFEST.sha256` lists the SHA-256 digest of every file in this directory except
itself. Verify with

```
python -c "import hashlib,pathlib;[print(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest(),p) for p in ['replay.py','README.md']]"
```

or with any standard sha256 checker that understands the `<digest>  <name>` format.
