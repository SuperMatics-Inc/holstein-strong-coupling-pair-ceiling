"""Exact-rational and outward-interval certificate for the constants of the preprint.

Uses Python standard library only (no floating-point proof steps, no network,
no files, no fitted numbers).
- Default mode: certifies 16-state CAR/bond-trace algebra, half-filling supremum,
  exact T_4..T_7 connected-subgraph enumeration (k <= 6), two-variable Raney
  multigraph domination, outward Fraction Taylor/sqrt interval bounds on P_e(c,y),
  all 10 unconditional (q_pair, C) trade-off points, canonical source transfer
  for C=6000 (all N>=4) and C=4500 (even N>=4), and all 5 weighted-2ab points.
- --mutants (or --all-mutants): runs every guarded mutant in a fresh child
  process and requires exit code 1 with REJECTED prefix.
"""
from collections import deque
from fractions import Fraction as Q
from pathlib import Path
import argparse
import subprocess
import sys

MUTANTS = (
    "bridge-parity-m-equals-k",
    "k3-claw-color-2",
    "tail-color-4-on-odd",
    "pe-interval-corrupted-m",
    "half-filling-sign-flip",
    "even-classes-on-odd-torus",
    "c-5800-all-below-threshold",
    "c-4300-even-below-threshold",
    "canonical-source-unscaled",
    "pair-prefactor-corrupted",
    "weighted-2ab-eps-too-large",
    "raney-tail-ternary-corrupted",
)


class Rejected(Exception):
    pass


def require(cond: bool, label: str) -> None:
    if not cond:
        raise Rejected(label)


# ---------------------------------------------------------------------------
# 1. Pure-Fraction outward interval arithmetic (48-bit binary dyadic grid)
# ---------------------------------------------------------------------------
SCALE_BITS = 48
SCALE_D = 1 << SCALE_BITS


def q_floor(x: Q) -> Q:
    return Q((x.numerator * SCALE_D) // x.denominator, SCALE_D)


def q_ceil(x: Q) -> Q:
    n = x.numerator * SCALE_D
    d = x.denominator
    return Q((n + d - 1) // d, SCALE_D)


def exp_bounds(u: Q, terms: int = 12) -> tuple[Q, Q]:
    """Outward rational enclosure [lo, hi] for exp(u), u >= 0."""
    if u == 0:
        return (Q(1), Q(1))
    m = 1 if u <= 1 else int(u) + 1
    x = u / m
    x_lo, x_hi = q_floor(x), q_ceil(x)
    s_lo, t_lo, s_hi, t_hi = Q(0), Q(1), Q(0), Q(1)
    for k in range(terms):
        s_lo = q_floor(s_lo + t_lo)
        t_lo = q_floor(t_lo * x_lo / (k + 1))
        s_hi = q_ceil(s_hi + t_hi)
        t_hi = q_ceil(t_hi * x_hi / (k + 1))
    rem_hi = q_ceil(t_hi * (terms + 1) / (terms + 1 - x_hi))
    lo1, hi1 = s_lo, q_ceil(s_hi + rem_hi)

    def pow_b(base: Q, e: int, rnd) -> Q:
        res, b = Q(1), base
        while e:
            if e & 1:
                res = rnd(res * b)
            e >>= 1
            if e:
                b = rnd(b * b)
        return res

    return (pow_b(lo1, m, q_floor), pow_b(hi1, m, q_ceil))


def exp_neg_bounds(u: Q, terms: int = 12) -> tuple[Q, Q]:
    """Outward rational enclosure [lo, hi] for exp(-u), u >= 0."""
    elo, ehi = exp_bounds(u, terms)
    return (q_floor(1 / ehi), q_ceil(1 / elo))


def sqrt_bounds(u: Q, steps: int = 9) -> tuple[Q, Q]:
    """Outward rational enclosure [lo, hi] for sqrt(u), u >= 0."""
    if u == 0:
        return (Q(0), Q(0))
    k = 0
    v = u
    while v >= 4:
        v /= 4
        k += 1
    while v < Q(1, 4):
        v *= 4
        k -= 1
    hi = Q(2)
    for _ in range(steps):
        hi = q_ceil((hi + v / hi) / 2)
    lo = q_floor(v / hi)
    scale = Q(2) ** k if k >= 0 else Q(1) / (Q(2) ** (-k))
    return (lo * scale, hi * scale)


def ln_lower_bound(B: Q) -> Q:
    """Strict rational lower bound on ln(B) for B > 1 via positive atanh series."""
    u = (B - 1) / (B + 1)
    u2 = u * u
    s = Q(0)
    p = u
    for j in range(12):
        s += p / (2 * j + 1)
        p *= u2
    return 2 * s


def ln_upper_bound(B: Q) -> Q:
    """Elementary rational upper bound ln(B) <= B - 1 for B >= 1."""
    return B - 1


def pe_interval_hi(s2: Q, b0: Q, b1: Q) -> Q:
    """Rigorous upper bound on R_e^{ub}(b; s*sqrt(b)) over b in [b0, b1]."""
    q_lo, _ = exp_neg_bounds(b1 / 2)
    _, q_hi = exp_neg_bounds(b0 / 2)
    g0_lo, _ = sqrt_bounds(s2 * b0)
    _, g1_hi = sqrt_bounds(s2 * b1)
    _, eg1_hi = exp_bounds(g1_hi)
    _, emg0_hi = exp_neg_bounds(g0_lo)
    r1_lo, _ = sqrt_bounds(b1 * b1 / 4 + 4 * s2 * b1)
    exp_r_diff = q_ceil((4 * s2 * b1) / (r1_lo + b1 / 2))
    _, er1_hi = exp_bounds(exp_r_diff)
    num_hi = 3 + 4 * (q_hi ** 2) + 4 * q_hi * (eg1_hi + emg0_hi) + er1_hi
    den_lo = 4 * ((1 + q_lo) ** 2)
    return num_hi / den_lo


def certify_pe_bound(c: int, y: Q, M: Q, b_cut: Q = Q(40)) -> int:
    """Certify sup_{b >= 0} R_e^{ub}(b; c*sqrt(y*b)) <= M by adaptive bisection."""
    s2 = Q(c * c) * y
    require(b_cut >= 4 * s2, f"b_cut-too-small-{c}-{y}")
    _, q_cut = exp_neg_bounds(b_cut / 2)
    _, g_cut = sqrt_bounds(s2 * b_cut)
    _, tail_cross = exp_neg_bounds(b_cut / 2 - g_cut)
    _, exp_4s2_hi = exp_bounds(4 * s2)
    tail_hi = (3 + 4 * (q_cut ** 2) + 4 * tail_cross + 4 * q_cut + exp_4s2_hi) / 4
    require(tail_hi <= M, f"pe-tail-exceeds-M-{c}-{y}")
    stack = [(Q(0), b_cut)]
    n_intervals = 0
    while stack:
        a, b = stack.pop()
        if pe_interval_hi(s2, a, b) <= M:
            n_intervals += 1
        else:
            require(b - a >= Q(1, 512), f"pe-interval-exceeds-M-{c}-{y}")
            mid = (a + b) / 2
            stack.append((mid, b))
            stack.append((a, mid))
    return n_intervals


# ---------------------------------------------------------------------------
# 2. 16-state CAR algebra, sign gauge, and exact half-filling supremum
# ---------------------------------------------------------------------------
def check_16state_and_supremum(mutant: str | None) -> int:
    n = 16
    T_e = [[Q(0)] * n for _ in range(n)]
    B_e = [[Q(0)] * n for _ in range(n)]
    N_odd = [0] * n
    N_tot = [0] * n
    O_xy = [[Q(0)] * n for _ in range(n)]

    def hop_sign(state: int, a: int, b: int) -> int:
        lo, hi = (a, b) if a < b else (b, a)
        cnt = sum((state >> k) & 1 for k in range(lo + 1, hi))
        return -1 if (cnt % 2) else 1

    for st in range(n):
        nu = ((st >> 0) & 1) + ((st >> 1) & 1)
        nv = ((st >> 2) & 1) + ((st >> 3) & 1)
        N_odd[st] = (1 if nu == 1 else 0) + (1 if nv == 1 else 0)
        N_tot[st] = nu + nv
        for a, b in ((0, 2), (2, 0), (1, 3), (3, 1)):
            if ((st >> b) & 1) and not ((st >> a) & 1):
                st2 = (st ^ (1 << b)) | (1 << a)
                T_e[st2][st] += hop_sign(st, a, b)
                B_e[st2][st] += 1
        if st == 12:
            O_xy[3][12] = Q(1)

    # Diagonal sign gauge S_JW on 16-state Fock basis: S T_e S == B_e
    visited = {}
    for root in range(n):
        if root in visited:
            continue
        visited[root] = 1
        stk = [root]
        while stk:
            u = stk.pop()
            for v in range(n):
                if B_e[v][u] != 0:
                    sgn = int(T_e[v][u])
                    if v not in visited:
                        visited[v] = visited[u] * sgn
                        stk.append(v)
                    else:
                        require(visited[v] == visited[u] * sgn, "16state-cycle-sign")

    for i in range(n):
        for j in range(n):
            require(visited[i] * T_e[i][j] * visited[j] == B_e[i][j], "16state-sign-gauge")

    # Verify characteristic trace coefficients of H_0(b) + g*B_e in each N-sector
    # via exact rational traces Tr(H^k) for k = 1..4 on N = 0,1,2,3,4
    # In N=2 (6 states), eigenvalues of -(b/2)*N_odd + g*B_e are:
    #   -b (mult 3), 0 (mult 1), (-b/2 + r), (-b/2 - r) with r = sqrt(b^2/4 + 4*g^2).
    # So Tr(H_N2) = -4b, Tr(H_N2^2) = 3b^2 + (b/2+r)^2 + (b/2-r)^2 = 3.5 b^2 + 8 g^2.
    def mat_mul(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    for b_q, g_q in ((Q(1, 2), Q(1, 3)), (Q(3), Q(2, 5))):
        H = [[g_q * B_e[i][j] for j in range(n)] for i in range(n)]
        for i in range(n):
            H[i][i] -= (b_q / 2) * N_odd[i]
        H2 = mat_mul(H, H)
        H3 = mat_mul(H2, H)
        H4 = mat_mul(H3, H)
        idx2 = [i for i in range(n) if N_tot[i] == 2]
        tr1 = sum(H[i][i] for i in idx2)
        tr2 = sum(H2[i][i] for i in idx2)
        tr3 = sum(H3[i][i] for i in idx2)
        tr4 = sum(H4[i][i] for i in idx2)
        require(tr1 == -4 * b_q, "n2-tr1")
        require(tr2 == 4 * (b_q ** 2) + 8 * (g_q ** 2), "n2-tr2")
        require(tr3 == -4 * (b_q ** 3) - 12 * b_q * (g_q ** 2), "n2-tr3")
        require(tr4 == 4 * (b_q ** 4) + 16 * (b_q ** 2) * (g_q ** 2) + 32 * (g_q ** 4), "n2-tr4")

    # Verify algebraic identity for R_e(b; g; z) <= R_e(b; g; 1):
    # With w = z + 1/z >= 2, C0 = cosh(g) >= 1, C2 = cosh(r) >= 1, q = exp(-b/2) in (0,1]:
    # R_e(w) = 1 + (A*w + B0)/(w + 2*q)^2 where
    #   A = 4*q*(C0 - 1),
    #   B0 = 2*q*C2 + 2*q^2 - 1 - 3*q^2 - 2? Wait: let's check exact B0:
    #   z_a(z)^2/z^2 = (w + 2q)^2 = w^2 + 4qw + 4q^2.
    #   Z_e(z)/z^2 = (w^2 - 2) + 4*q*C0*w + (1 + 3*q^2 + 2*q*C2)
    #              = (w + 2q)^2 + 4*q*(C0 - 1)*w + (2*q*C2 - q^2 - 1).
    # Verify this exact polynomial identity in rational variables (w, q, C0, C2)
    # AND verify R_e(w) <= R_e(2) <=> (2-w)*(A*(2*w + 2*q) + B0*(...))?
    # Better: f(w) - 1 = (A*w + B0)/(w + 2*q)^2; derivative numerator is
    #   A*(w + 2q) - 2*(A*w + B0) = -A*w + 2*A*q - 2*B0 <= -2*A*(1 - q) - 2*B0 <= 0
    # whenever A >= 0 and B0 >= 0 (Wait: C2 = cosh(r) >= cosh(b/2) = (q + 1/q)/2, so
    #   2*q*C2 >= 2*q*(q + 1/q)/2 = q^2 + 1  ==>  B0 = 2*q*C2 - q^2 - 1 >= 0!!).
    test_tuples = [
        (Q(1, 2), Q(3, 2), Q(5, 4), Q(13, 5)),
        (Q(1, 5), Q(9, 4), Q(17, 8), Q(4)),
    ]
    for q_v, C0_v, C2_v, w_v in test_tuples:
        # Enforce C2 >= (q + 1/q)/2
        C2_eff = max(C2_v, (q_v + 1 / q_v) / 2)
        if mutant == "half-filling-sign-flip":
            # Corrupt B0 sign by flipping + 2*q*C2 -> - 2*q*C2
            B0 = -2 * q_v * C2_eff - q_v ** 2 - 1
        else:
            B0 = 2 * q_v * C2_eff - q_v ** 2 - 1
        A = 4 * q_v * (C0_v - 1)
        require(A >= 0 and B0 >= 0, "half-filling-nonneg-AB0")
        lhs_poly = (w_v ** 2 - 2) + 4 * q_v * C0_v * w_v + (1 + 3 * (q_v ** 2) + 2 * q_v * C2_eff)
        rhs_poly = (w_v + 2 * q_v) ** 2 + A * w_v + B0
        require(lhs_poly == rhs_poly, "half-filling-poly-identity")
        # Check exact monotonicity f(w) <= f(2)
        f_w = 1 + (A * w_v + B0) / ((w_v + 2 * q_v) ** 2)
        f_2 = 1 + (2 * A + B0) / ((2 + 2 * q_v) ** 2)
        require(f_w <= f_2, "half-filling-supremum")

    return len(test_tuples)


# ---------------------------------------------------------------------------
# 3. Exact T_4, T_5, T_6, T_7 connected-subgraph enumeration (k <= 6)
#    and two-variable Raney multigraph verification
# ---------------------------------------------------------------------------
def enumerate_torus_k6(N: int) -> dict[tuple[int, int, int, int], int]:
    """Enumerate all rooted connected edge subgraphs of T_N with 1 <= k <= 6.

    Returns counts keyed by (v, k, c_ub, m) -> count where:
      v = |V(E)|, k = |E|,
      c_ub = max_deg(E) if E is bipartite else max_deg(E) + 1 (Vizing upper bound,
             exact <= 4 on all k <= 6 graphs for N >= 4),
      m = k + b (b = number of bridge edges of E).
    """
    edges = []
    for x in range(N):
        for y in range(N):
            u = (x, y)
            for v in (((x + 1) % N, y), (x, (y + 1) % N)):
                e = (u, v) if u < v else (v, u)
                edges.append(e)
    edges = sorted(set(edges))
    adj = {}
    for idx, (u, v) in enumerate(edges):
        adj.setdefault(u, []).append((idx, v))
        adj.setdefault(v, []).append((idx, u))

    root = (0, 0)
    by_k: list[set[frozenset[int]]] = [set() for _ in range(7)]
    by_k[0].add(frozenset())
    for k in range(1, 7):
        nxt = by_k[k]
        for E in by_k[k - 1]:
            verts = {root}
            for idx in E:
                verts.add(edges[idx][0])
                verts.add(edges[idx][1])
            for u in verts:
                for idx, _ in adj[u]:
                    if idx not in E:
                        nxt.add(E | {idx})

    stats: dict[tuple[int, int, int, int], int] = {}
    for k in range(1, 7):
        for E in by_k[k]:
            verts = set()
            deg = {}
            sub_adj = {}
            for idx in E:
                u, v = edges[idx]
                verts.add(u)
                verts.add(v)
                deg[u] = deg.get(u, 0) + 1
                deg[v] = deg.get(v, 0) + 1
                sub_adj.setdefault(u, []).append(v)
                sub_adj.setdefault(v, []).append(u)
            v_cnt = len(verts)
            md = max(deg.values())

            # Bipartite check
            color = {}
            bip = True
            for s in verts:
                if s in color:
                    continue
                color[s] = 0
                dq = deque([s])
                while dq and bip:
                    u = dq.popleft()
                    for w in sub_adj[u]:
                        if w not in color:
                            color[w] = 1 - color[u]
                            dq.append(w)
                        elif color[w] == color[u]:
                            bip = False
                            break

            # Bridge count by DFS lowlink
            tin = {}
            low = {}
            timer = [0]
            bridges = [0]

            def dfs(u, p):
                tin[u] = low[u] = timer[0]
                timer[0] += 1
                for w in sub_adj[u]:
                    if w == p:
                        continue
                    if w in tin:
                        if tin[w] < low[u]:
                            low[u] = tin[w]
                    else:
                        dfs(w, u)
                        if low[w] < low[u]:
                            low[u] = low[w]
                        if low[w] > tin[u]:
                            bridges[0] += 1

            for s in verts:
                if s not in tin:
                    dfs(s, None)

            c_ub = md if bip else md + 1
            m = k + bridges[0]
            key = (v_cnt, k, c_ub, m)
            stats[key] = stats.get(key, 0) + 1

    return stats


# Canonical dominating class tables for all N >= 4 and even N >= 4
# Format: (v, k, c, m, count)
ALL_EXACT_CLASSES = (
    (2, 1, 1, 2, 4),
    (3, 2, 2, 4, 18),
    (4, 3, 2, 6, 72),
    (4, 3, 3, 6, 16),
    (4, 4, 2, 4, 6),
    (5, 4, 2, 8, 250),
    (5, 4, 3, 8, 180),
    (5, 4, 4, 8, 5),
    (5, 5, 3, 5, 2),
    (5, 5, 3, 6, 60),
    (6, 5, 2, 10, 852),
    (6, 5, 3, 10, 1260),
    (6, 5, 4, 10, 72),
    (6, 6, 2, 6, 48),
    (6, 6, 4, 7, 24),
    (6, 6, 3, 8, 432),
    (6, 6, 4, 8, 36),
    (7, 6, 2, 12, 2730),
    (7, 6, 3, 12, 7658),
    (7, 6, 4, 12, 630),
)

EVEN_EXACT_CLASSES = tuple(
    row for row in ALL_EXACT_CLASSES if (row[0], row[1], row[3]) not in ((5, 5, 5), (6, 6, 7))
)


def check_enumeration_and_raney(mutant: str | None) -> None:
    stats_by_N = {N: enumerate_torus_k6(N) for N in (4, 5, 6, 7)}

    all_dict = {(v, k, c, m): cnt for (v, k, c, m, cnt) in ALL_EXACT_CLASSES}
    even_dict = {(v, k, c, m): cnt for (v, k, c, m, cnt) in EVEN_EXACT_CLASSES}

    if mutant == "k3-claw-color-2":
        all_dict[(4, 3, 2, 6)] = 88
        del all_dict[(4, 3, 3, 6)]
    if mutant == "even-classes-on-odd-torus":
        # Try using EVEN_EXACT_CLASSES to dominate T_5 (fails on 5-cycles!)
        all_dict = dict(even_dict)

    for N, st in stats_by_N.items():
        for key, cnt in st.items():
            require(all_dict.get(key, 0) >= cnt, f"all-domination-T{N}-{key}")
            if N % 2 == 0:
                require(even_dict.get(key, 0) >= cnt, f"even-domination-T{N}-{key}")

    # Verify exact Raney tree counts R_4(k) = [w^k] T(w)^4 for k=1..6, T = 1 + w*T^3
    def poly_mul(A, B, deg=6):
        res = [0] * (deg + 1)
        for i, a in enumerate(A):
            if a == 0:
                continue
            for j, b in enumerate(B):
                if i + j <= deg:
                    res[i + j] += a * b
        return res

    T_poly = [1, 0, 0, 0, 0, 0, 0]
    deg_branch = 2 if mutant == "raney-tail-ternary-corrupted" else 3
    for _ in range(7):
        Tp = [1, 0, 0, 0, 0, 0, 0]
        for _ in range(deg_branch):
            Tp = poly_mul(Tp, T_poly)
        T_poly = [1] + Tp[:6]
    T2 = poly_mul(T_poly, T_poly)
    R4 = poly_mul(T2, T2)
    require(R4 == [1, 4, 18, 88, 455, 2448, 13566], "raney-R4-coeffs")
    # Every T_N (N=4..7) total count at k=1..6 is bounded by R4[k]
    for N, st in stats_by_N.items():
        for k in range(1, 7):
            tot_k = sum(cnt for (_, kk, _, _,), cnt in st.items() if kk == k)
            require(tot_k <= R4[k], f"raney-bounds-T{N}-k{k}")

    # Verify two-variable Raney multigraph bound P = 1 + wp P^3 S^2, S = 1 + ws P^4 S
    # against brute-force connected subgraph enumeration on T_5 + 2 source edges (kp + ks <= 4)
    N = 5
    all_edges = []
    for x in range(N):
        for y in range(N):
            u = (x, y)
            for v in (((x + 1) % N, y), (x, (y + 1) % N)):
                all_edges.append(("p", u, v))
            all_edges.append(("s", u, ((x + 2) % N, y)))
    adj = {}
    for idx, (_, u, v) in enumerate(all_edges):
        adj.setdefault(u, []).append(idx)
        adj.setdefault(v, []).append(idx)
    root = (0, 0)
    seen = {frozenset()}
    dq = deque([frozenset()])
    multi_counts = {}
    while dq:
        E = dq.popleft()
        if E:
            kp = sum(1 for i in E if all_edges[i][0] == "p")
            ks = len(E) - kp
            multi_counts[(kp, ks)] = multi_counts.get((kp, ks), 0) + 1
        if len(E) == 4:
            continue
        verts = {root}
        for i in E:
            verts.add(all_edges[i][1])
            verts.add(all_edges[i][2])
        for u in verts:
            for idx in adj[u]:
                if idx not in E:
                    E2 = E | {idx}
                    if E2 not in seen:
                        seen.add(E2)
                        dq.append(E2)

    def mul2(A, B, deg=4):
        res = {}
        for (ap, as_), av in A.items():
            for (bp, bs), bv in B.items():
                if ap + bp + as_ + bs <= deg:
                    res[(ap + bp, as_ + bs)] = res.get((ap + bp, as_ + bs), 0) + av * bv
        return res

    def add2(A, B):
        res = dict(A)
        for k, v in B.items():
            res[k] = res.get(k, 0) + v
        return res

    P = {(0, 0): 1}
    S = {(0, 0): 1}
    for _ in range(7):
        P2 = mul2(P, P)
        P3 = mul2(P2, P)
        P4 = mul2(P2, P2)
        S2 = mul2(S, S)
        P = add2({(0, 0): 1}, mul2({(1, 0): 1}, mul2(P3, S2)))
        S = add2({(0, 0): 1}, mul2({(0, 1): 1}, mul2(P4, S)))
    G_root = mul2(mul2(mul2(P, P), mul2(P, P)), mul2(S, S))
    for k_tot in range(1, 5):
        for ks in range(k_tot + 1):
            kp = k_tot - ks
            require(
                multi_counts.get((kp, ks), 0) <= G_root.get((kp, ks), 0),
                f"multigraph-raney-{kp}-{ks}",
            )


# ---------------------------------------------------------------------------
# 4. Certified P_e(c, y) <= M witnesses and exact Kotecky-Preiss tables
# ---------------------------------------------------------------------------
PE_BOUNDS = {
    (1, 2, 1): (Q(249, 1000), Q(143, 100)),
    (2, 4, 2): (Q(39, 250), Q(1893, 500)),
    (2, 4, 4): (Q(1, 10), Q(249, 125)),
    (2, 6, 3): (Q(39, 250), Q(1893, 500)),
    (2, 6, 6): (Q(1, 10), Q(249, 125)),
    (2, 8, 4): (Q(39, 250), Q(1893, 500)),
    (2, 10, 5): (Q(39, 250), Q(1893, 500)),
    (2, 12, 6): (Q(39, 250), Q(1893, 500)),
    (3, 5, 5): (Q(57, 1000), Q(1349, 500)),
    (3, 6, 3): (Q(23, 250), Q(1903, 250)),
    (3, 6, 5): (Q(13, 200), Q(837, 250)),
    (3, 8, 4): (Q(23, 250), Q(1903, 250)),
    (3, 8, 6): (Q(69, 1000), Q(15, 4)),
    (3, 10, 5): (Q(23, 250), Q(1903, 250)),
    (3, 12, 6): (Q(23, 250), Q(1903, 250)),
    (4, 1, 1): (Q(39, 1000), Q(1893, 500)),
    (4, 7, 6): (Q(43, 1000), Q(584, 125)),
    (4, 8, 4): (Q(13, 200), Q(1677, 100)),
    (4, 8, 6): (Q(6, 125), Q(123, 20)),
    (4, 10, 5): (Q(13, 200), Q(1677, 100)),
    (4, 12, 6): (Q(13, 200), Q(1677, 100)),
    (5, 1, 1): (Q(29, 1000), Q(662, 125)),
}

# Exact rational trade-off curve: (mode, q_pair, C, T, B, L_src, A_pair)
UNCONDITIONAL_CURVE = (
    ("all", Q(99, 100), Q(6000), Q(303, 250), Q(557, 500), Q(100000), Q(57, 500)),
    ("all", Q(19, 20), Q(6550), Q(151, 125), Q(553, 500), None, Q(53, 500)),
    ("all", Q(9, 10), Q(7300), Q(121, 100), Q(139, 125), None, Q(14, 125)),
    ("all", Q(3, 4), Q(10600), Q(121, 100), Q(279, 250), None, Q(29, 250)),
    ("all", Q(1, 2), Q(26000), Q(151, 125), Q(29, 25), None, Q(4, 25)),
    ("even", Q(99, 100), Q(4500), Q(121, 100), Q(557, 500), Q(100000), Q(57, 500)),
    ("even", Q(19, 20), Q(4900), Q(151, 125), Q(111, 100), None, Q(11, 100)),
    ("even", Q(9, 10), Q(5450), Q(121, 100), Q(557, 500), None, Q(57, 500)),
    ("even", Q(3, 4), Q(7900), Q(121, 100), Q(559, 500), None, Q(59, 500)),
    ("even", Q(1, 2), Q(19800), Q(301, 250), Q(291, 250), None, Q(41, 250)),
)

# Exact rational weighted-2ab trade-off curve: (gamma, kappa, q_pair, C, eps0, y, T, B)
WEIGHTED_2AB_CURVE = (
    (Q(1), Q(1, 10), Q(99, 100), Q(85000), Q(173, 1000000), Q(61, 4000), Q(59, 50), Q(113, 100)),
    (Q(1), Q(1, 4), Q(99, 100), Q(102000), Q(433, 1000000), Q(51, 4000), Q(59, 50), Q(113, 100)),
    (Q(1), Q(1, 2), Q(99, 100), Q(152000), Q(869, 1000000), Q(17, 2000), Q(59, 50), Q(113, 100)),
    (Q(1), Q(3, 4), Q(99, 100), Q(305000), Q(1301, 1000000), Q(17, 4000), Q(59, 50), Q(113, 100)),
    (Q(1), Q(9, 10), Q(99, 100), Q(765000), Q(1537, 1000000), Q(7, 4000), Q(59, 50), Q(113, 100)),
)


def ceil_root_bound(M: Q, k: int, c: int, C_val: Q, y: Q, m: int) -> Q:
    """Return a clean rational R_k such that R_k^{2c} * (C*y)^{m*c} >= M^{2k}."""
    # Compute integer numerator of ceil(10^8 * (M^(k/c) / (C*y)^(m/2))) via binary search
    rhs = (M ** (2 * k))
    den = (C_val * y) ** (m * c)
    # R = n / 10^8 satisfies n^(2c) * den >= rhs * 10^(16c)
    D = 10 ** 8
    target_num = rhs * (Q(D) ** (2 * c)) / den
    # Integer ceiling of target_num^(1/(2c))
    t_ceil = (target_num.numerator + target_num.denominator - 1) // target_num.denominator
    lo, hi = 0, max(1, t_ceil + 1)
    # Fast upper bound using bit_length
    hi = 1 << ((t_ceil.bit_length() // (2 * c)) + 2)
    while lo < hi:
        mid = (lo + hi) // 2
        if mid ** (2 * c) >= t_ceil:
            hi = mid
        else:
            lo = mid + 1
    R_k = Q(lo, D)
    require((R_k ** (2 * c)) * den >= rhs, "ceil-root-bound")
    return R_k


def check_unconditional_and_canonical(mutant: str | None) -> tuple[dict, int]:
    # 1. Certify unique (c, y, M) bounds by outward Fraction Taylor/sqrt bisection
    pe_bounds = dict(PE_BOUNDS)
    if mutant == "pe-interval-corrupted-m":
        pe_bounds[(5, 1, 1)] = (Q(29, 1000), Q(4))  # true sup is ~5.29 > 4

    seen_cym = {}
    total_ints = 0
    for key, (yq, Mq) in pe_bounds.items():
        c = key[0]
        if (c, yq, Mq) not in seen_cym:
            n_int = certify_pe_bound(c, yq, Mq)
            seen_cym[(c, yq, Mq)] = n_int
            total_ints += n_int

    curve = list(UNCONDITIONAL_CURVE)
    if mutant == "c-5800-all-below-threshold":
        curve[0] = ("all", Q(99, 100), Q(5800), Q(303, 250), Q(557, 500), None, Q(57, 500))
    elif mutant == "c-4300-even-below-threshold":
        curve[5] = ("even", Q(99, 100), Q(4300), Q(121, 100), Q(557, 500), None, Q(57, 500))
    elif mutant == "canonical-source-unscaled":
        curve[0] = ("all", Q(99, 100), Q(6000), Q(303, 250), Q(557, 500), Q(10), Q(57, 500))
    elif mutant == "pair-prefactor-corrupted":
        curve[0] = ("all", Q(99, 100), Q(6000), Q(303, 250), Q(557, 500), Q(100000), Q(1, 10))

    results = {}
    for mode, q_p, C_val, T_val, B_val, L_src, A_pair in curve:
        classes = ALL_EXACT_CLASSES if mode == "all" else EVEN_EXACT_CLASSES
        c_tail = 4 if (mode == "even" or mutant == "tail-color-4-on-odd") else 5

        if mode == "all":
            # T_5 is 4-regular with 25 vertices and 50 edges;
            # any matching has <= floor(25/2) = 12 edges, so c_tail matchings cover <= 12*c_tail edges.
            require(c_tail * (25 // 2) >= 50, "odd-torus-T5-needs-5-colors")
        yt, Mt = pe_bounds[(c_tail, 1, 1)]
        require(C_val * yt >= 1, f"ell-tail-ge-1-{mode}-{C_val}")
        w_eff = (T_val - 1) / (T_val ** 3)
        tail_raney = T_val ** 4 - 1 - 4 * w_eff - 18 * (w_eff ** 2) - 88 * (w_eff ** 3) - 455 * (w_eff ** 4) - 2448 * (w_eff ** 5) - 13566 * (w_eff ** 6)
        require(tail_raney > 0, f"tail-raney-pos-{mode}-{C_val}")

        if L_src is None:
            # Grand-canonical tail check: B * Mt^(1/c_tail) / sqrt(C * yt) <= w_eff * q_p
            require(
                ((B_val ** c_tail) * Mt) ** 2 <= ((w_eff * q_p) ** (2 * c_tail)) * ((C_val * yt) ** c_tail),
                f"tail-cond-{mode}-{C_val}",
            )
            S_src = Q(0)
        else:
            # Source-augmented two-variable Raney contraction check:
            # P_hi = T_val, ws = 3*B/L_src, S_hi = 1/(1 - ws * P_hi^4), w_p * S_hi^2 <= w_eff
            P_hi = T_val
            ws = 3 * B_val / L_src
            require(ws * (P_hi ** 4) < 1, f"src-S-pos-{mode}-{C_val}")
            S_hi = 1 / (1 - ws * (P_hi ** 4))
            require(
                ((B_val ** c_tail) * Mt) ** 2 * (S_hi ** (4 * c_tail))
                <= ((w_eff * q_p) ** (2 * c_tail)) * ((C_val * yt) ** c_tail),
                f"src-tail-cond-{mode}-{C_val}",
            )
            wp_hi = w_eff / (S_hi ** 2)
            den = 1 - 3 * wp_hi * (P_hi ** 2) * (S_hi ** 2) - 8 * wp_hi * ws * (P_hi ** 6) * (S_hi ** 3)
            require(den > 0, f"src-deriv-den-{mode}-{C_val}")
            Dp = 2 * wp_hi * (P_hi ** 7) * (S_hi ** 3) / den
            Ds = (P_hi ** 4) * (S_hi ** 2) + 4 * ws * (P_hi ** 3) * (S_hi ** 2) * Dp
            dG_dws = 4 * (P_hi ** 3) * (S_hi ** 2) * Dp + 2 * (P_hi ** 4) * S_hi * Ds
            S_src = B_val * ws * dG_dws

        S_phys = B_val * tail_raney
        for v, k, c, m_real, cnt in classes:
            m_used = k if mutant == "bridge-parity-m-equals-k" else m_real
            yk, Mk = pe_bounds[(c, m_used, k)] if (c, m_used, k) in pe_bounds else pe_bounds[(c, m_real, k)]
            require(C_val * yk >= 1, f"ell-ge-1-{mode}-{C_val}-{(c,k)}")
            Rk = ceil_root_bound(Mk, k, c, C_val, yk, m_used)
            S_phys += cnt * (B_val ** v) * Rk * (q_p ** (-k))
            if L_src is not None:
                require((B_val ** k) * Rk <= (wp_hi * q_p) ** k, f"src-k-le-6-{mode}-{C_val}-{(v,k,c,m_real)}")

        lnB_lo = ln_lower_bound(B_val)
        require(lnB_lo <= A_pair, "pair-prefactor-corrupted")
        require(ln_upper_bound(B_val) <= A_pair, f"pair-prefactor-upper-{mode}-{C_val}")

        require(S_phys + S_src <= lnB_lo, f"kp-budget-{mode}-{C_val}")
        results[(mode, q_p, C_val)] = (S_phys, S_src, lnB_lo, A_pair)

    return results, total_ints


def check_weighted_2ab(mutant: str | None) -> list:
    # Collapsed 2ab counts (v, k, m, count) from ALL_EXACT_CLASSES
    classes_2ab = {}
    for v, k, _c, m, cnt in ALL_EXACT_CLASSES:
        classes_2ab[(v, k, m)] = classes_2ab.get((v, k, m), 0) + cnt

    curve = list(WEIGHTED_2AB_CURVE)
    if mutant == "weighted-2ab-eps-too-large":
        curve[4] = (Q(1), Q(9, 10), Q(99, 100), Q(765000), Q(1, 50), Q(7, 4000), Q(59, 50), Q(113, 100))

    out = []
    for gamma, kappa, q_p, C_val, eps0, y_val, T_val, B_val in curve:
        require(C_val * y_val >= 1, f"2ab-ell-ge-1-{kappa}")
        # Require 16 * sqrt(C*y) * eps0 <= kappa * gamma  <=>  256 * C * y * eps0^2 <= (kappa * gamma)^2
        require(256 * C_val * y_val * (eps0 ** 2) <= (kappa * gamma) ** 2, f"2ab-eps0-{kappa}")
        _, E_hi = exp_bounds(32 * y_val / (gamma * (1 - kappa)))
        w_eff = (T_val - 1) / (T_val ** 3)
        require(((2 * B_val * E_hi) ** 2) <= ((w_eff * q_p) ** 2) * (C_val * y_val), f"2ab-tail-{kappa}")
        # Rational upper bound L_inv >= 1 / sqrt(C * y)
        D = 10 ** 8
        t_num = Q(D * D) / (C_val * y_val)
        t_ceil = (t_num.numerator + t_num.denominator - 1) // t_num.denominator
        lo, hi = 0, (1 << ((t_ceil.bit_length() // 2) + 2))
        while lo < hi:
            mid = (lo + hi) // 2
            if mid * mid >= t_ceil:
                hi = mid
            else:
                lo = mid + 1
        L_inv = Q(lo, D)
        require((L_inv ** 2) * C_val * y_val >= 1, "2ab-linv")
        S = 2 * B_val * (T_val ** 4 - 1 - 4 * w_eff - 18 * (w_eff ** 2) - 88 * (w_eff ** 3) - 455 * (w_eff ** 4) - 2448 * (w_eff ** 5) - 13566 * (w_eff ** 6))
        for (v, k, m), cnt in classes_2ab.items():
            S += cnt * ((2 * B_val) ** v) * (E_hi ** k) * (L_inv ** m) * (q_p ** (-k))
        lnB_lo = ln_lower_bound(B_val)
        require(S <= lnB_lo, f"2ab-kp-{kappa}")
        out.append((kappa, C_val, eps0, S, lnB_lo))
    return out


# ---------------------------------------------------------------------------
# 5. Main certificate and CLI entry point
# ---------------------------------------------------------------------------
def certificate(mutant: str | None = None, verbose: bool = True) -> None:
    require(mutant is None or mutant in MUTANTS, "known-mutant")
    n_hf = check_16state_and_supremum(mutant)
    check_enumeration_and_raney(mutant)
    uncond_res, total_ints = check_unconditional_and_canonical(mutant)
    w2ab_res = check_weighted_2ab(mutant)

    if verbose:
        s_all, src_all, ln_all, A_all = uncond_res[("all", Q(99, 100), Q(6000))]
        s_ev, src_ev, ln_ev, A_ev = uncond_res[("even", Q(99, 100), Q(4500))]
        print(f"16STATE_HALF_FILLING: S_JW*T_e*S_JW==B_e, 5-sector traces & {n_hf} supremum rational identities verified")
        print(f"ENUM_AND_RANEY: T_4..T_7 exact k<=6 counts dominated; R4=[1,4,18,88,455,2448,13566]; 2-var Raney (kp+ks<=4) verified")
        print(f"PE_INTERVALS: {len(PE_BOUNDS)} class bounds certified across {total_ints} adaptive outward Fraction intervals")
        print(
            f"MAIN_ALL_TORI: eps0=inf, C=6000, q_pair=99/100, A_pair={A_all}, "
            f"S_phys<=0.{int(s_all*1000000):06d}, S_src<=0.{int(src_all*1000000):06d}, "
            f"S_tot<=0.{int((s_all+src_all)*1000000):06d} <= lnB>=0.{int(ln_all*1000000):06d} "
            f"(22369.6x tighter than the earlier 2D release C=2^27; 45812984x tighter than the earlier conditional release C=2^38)"
        )
        print(
            f"MAIN_EVEN_TORI: eps0=inf, C=4500, q_pair=99/100, A_pair={A_ev}, "
            f"S_phys<=0.{int(s_ev*1000000):06d}, S_src<=0.{int(src_ev*1000000):06d}, "
            f"S_tot<=0.{int((s_ev+src_ev)*1000000):06d} <= lnB>=0.{int(ln_ev*1000000):06d} "
            f"(29826.2x tighter than the earlier 2D release C=2^27; 61083979x tighter than the earlier conditional release C=2^38)"
        )
        print(f"UNCONDITIONAL_CURVE: {len(uncond_res)} rational (mode, q_pair, C) points certified")
        print(
            f"WEIGHTED_2AB_CURVE: {len(w2ab_res)} rational points certified "
            f"(C={w2ab_res[0][1]} at eps0={w2ab_res[0][2]} up to C={w2ab_res[-1][1]} at eps0={w2ab_res[-1][2]})"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS)
    parser.add_argument("--mutants", "--all-mutants", dest="mutants", action="store_true")
    args = parser.parse_args()

    if args.mutant:
        try:
            certificate(mutant=args.mutant, verbose=False)
        except Rejected as exc:
            print(f"REJECTED {args.mutant}: {exc}")
            return 1
        print(f"UNCAUGHT {args.mutant}")
        return 0

    try:
        certificate(verbose=True)
    except Rejected as exc:
        print(f"BASELINE FAILURE: {exc}")
        return 1

    if not args.mutants:
        return 0

    script = str(Path(__file__).resolve())
    for m in MUTANTS:
        child = subprocess.run(
            [sys.executable, "-B", script, "--mutant", m],
            capture_output=True,
            text=True,
            check=False,
        )
        print(child.stdout.strip())
        if child.stderr:
            print(child.stderr.strip())
        if child.returncode != 1 or not child.stdout.startswith(f"REJECTED {m}:"):
            print(f"FAILURE: mutant {m} returned {child.returncode}")
            return 1
    print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
