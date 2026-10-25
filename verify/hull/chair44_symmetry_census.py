#!/usr/bin/env python3
"""Which symmetry groups occur for Chair44 tilings? (Tsiokos, arXiv:2609.19214, Open Problem 1.)

Uses the 168-label cell substitution of the six-birds-tiles repository (via its own module
verify/substitution_modular_coincidence.py). Usage: python3 verify/hull/chair44_symmetry_census.py

A tiling with a symmetry group K fixing a vertex c of every supertile lattice is determined by the
2x2x2 label blocks w_n around c at each level n, with w_n = F(w_{n+1}) where F acts cellwise:
F(w)_i = phi(w_i, residue_i). Since F is R-equivariant, Stab(w_{n+1}) <= Stab(w_n) and the blocks
eventually cycle, so the groups that occur at such centres are exactly the stabilizers of
F-periodic blocks. Periodic blocks of the legal language give tilings in the substitution hull,
hence tilings by Q.
"""
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]  # repository root
sys.path.insert(0, str(REPO / "verify"))
import substitution_modular_coincidence as m  # noqa: E402

BITS = tuple(itertools.product((0, 1), repeat=3))
POS = tuple(itertools.product((-1, 0), repeat=3))
cert = json.loads(m.CERTIFICATE.read_text())
frames = tuple(m.frame_from_json(r) for r in cert["orientation_group"])
children = tuple((m.frame_from_json(r[0]), tuple(r[1])) for r in cert["children"])
labels = [(f, c) for f in frames for c in m.CHAIR_CELLS]
index = {l: i for i, l in enumerate(labels)}

phi = {}
for frame in frames:
    owners = {}
    for cf, co in m.child_poses(frame, children):
        for cell in m.CHAIR_CELLS:
            owners.setdefault(m.cell_lower(cf, co, cell), []).append((cf, cell))
    for cell in m.CHAIR_CELLS:
        src = m.cell_lower(frame, m.ORIGIN, cell)
        for d in BITS:
            (tl,) = owners[m.add(m.scale(2, src), d)]
            phi[index[(frame, cell)], d] = index[tl]

windows = [tuple((BITS.index(tuple((o[k] + b[k]) // 2 for k in range(3))),
                  tuple((o[k] + b[k]) % 2 for k in range(3))) for b in BITS)
           for o in itertools.product(range(3), repeat=3)]
language = {tuple(phi[j, d] for d in BITS) for j in range(168)}
frontier = set(language)
while frontier:
    new = {tuple(phi[blk[i], d] for i, d in win) for blk in frontier for win in windows} - language
    language |= new
    frontier = new
print("legal 2x2x2 blocks:", len(language))

residue = tuple(tuple(x % 2 for x in p) for p in POS)
# block indexed by BITS order in the language; seeds use POS order. Translate: language blocks are on
# {0,1}^3; the window with o = (1,1,1) of sigma(block) is centred on a vertex. Work with POS-blocks:
# the POS-blocks occurring in the language are the language blocks placed on {-1,0}^3 (same shape).
def F(w):
    return tuple(phi[w[i], residue[i]] for i in range(8))


def lower(R, x):
    return m.cell_lower(R, m.ORIGIN, x)


def act(R, w):
    conf = {lower(R, x): index[(m.compose(R, labels[l][0]), labels[l][1])] for x, l in zip(POS, w)}
    return tuple(conf[p] for p in POS)


assert all(F(w) in language for w in language), "F must preserve the language"

# periodic points of F on the language
periodic = set()
for w in language:
    seen = []
    x = w
    while x not in seen:
        seen.append(x)
        x = F(x)
    periodic.update(seen[seen.index(x):])
periods = Counter()
for w in periodic:
    p, x = 1, F(w)
    while x != w:
        x, p = F(x), p + 1
    periods[p] += 1
print("F-periodic legal blocks:", len(periodic), "by period:", dict(sorted(periods.items())))


def mat(R):
    return tuple(tuple(m.act(R, e)) for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))


def kind(R):
    """Classify a proper cube rotation: 'id', 'C2ax', 'C4', 'C2diag', 'C3'."""
    M = [list(c) for c in mat(R)]  # columns
    tr = M[0][0] + M[1][1] + M[2][2]
    def mul(A, B):
        return [[sum(A[k][i] * B[j][k] for k in range(3)) for i in range(3)] for j in range(3)]
    if tr == 3:
        return "id"
    if tr == 0:
        return "C3"
    if tr == 1:
        return "C4"
    # tr == -1: 180-degree rotation; axis along a coordinate iff a diagonal entry is +1
    return "C2ax" if any(M[i][i] == 1 for i in range(3)) else "C2diag"


def group_type(G):
    c = Counter(kind(R) for R in G)
    n = len(G)
    if n == 1:
        return "1"
    if n == 2:
        return "C2(axis)" if c["C2ax"] else "C2(diag)"
    if n == 4:
        if c["C4"]:
            return "C4"
        return "V4(axes)" if c["C2ax"] == 3 else "V4(axis+diags)"
    if n == 8:
        return "D4"
    return f"order {n} {dict(c)}"


stabs = {w: [R for R in frames if act(R, w) == w] for w in periodic}
hist = Counter((group_type(stabs[w]), ) for w in periodic)
print("symmetry groups of vertex-centred periodic tilings (block counts):")
for k, v in sorted(hist.items()):
    print("  ", k[0], v)
# also stabilizers of all legal blocks, for reference
allhist = Counter(group_type([R for R in frames if act(R, w) == w]) for w in language)
print("stabilizers over all legal blocks:", dict(allhist))
examples = {}
for w in sorted(periodic):
    examples.setdefault(group_type(stabs[w]), w)
print("one example block per group:", examples)

# ---------------------------------------------------------------- axial groups C2(axis), C2(diag), C4
# A cyclic symmetry group K with axis A: A passes through vertices of every supertile lattice; the
# blocks w_n at the vertex v_n of level n nearest a base point satisfy w_n = W_j(w_{n+1}), where W_j is
# the 2x2x2 window of sigma(w_{n+1}) centred at j*e (j in {-1,0,1}, e the axis direction). K occurs as
# the exact symmetry group of a hull tiling iff the graph of legal blocks with stabilizer exactly K
# (edges W_j) has a strongly connected component carrying two different shifts j (a constant shift
# sequence produces a vertex of every level on A, hence a vertex-centred tiling, whose group is 1 or D4).
LBITS = BITS


def window(w, o):
    return tuple(phi[w[BITS.index(tuple((o[k] + b[k]) // 2 for k in range(3)))],
                     tuple((o[k] + b[k]) % 2 for k in range(3))] for b in BITS)


assert all(window(w, (1, 1, 1)) == F(w) for w in language)
ident = [R for R in frames if all(m.act(R, e) == e for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))][0]


def axis_of(R):
    return [e for e in itertools.product((-1, 0, 1), repeat=3) if any(e) and m.act(R, e) == e]


def scc(nodes, edges):
    idx, low, st, on, out, c = {}, {}, [], set(), [], [0]
    sys.setrecursionlimit(100000)

    def visit(v):
        idx[v] = low[v] = c[0]; c[0] += 1; st.append(v); on.add(v)
        for u, _ in edges.get(v, ()):
            if u not in idx:
                visit(u); low[v] = min(low[v], low[u])
            elif u in on:
                low[v] = min(low[v], idx[u])
        if low[v] == idx[v]:
            comp = set()
            while True:
                u = st.pop(); on.discard(u); comp.add(u)
                if u == v:
                    break
            out.append(comp)
    for v in nodes:
        if v not in idx:
            visit(v)
    return out


stab_all = {w: frozenset(R for R in frames if act(R, w) == w) for w in language}
# For a cyclic group K about an axis with direction e, only the rotations about e commute with the axial
# windows, so the blocks to follow are those whose stabiliser meets the rotations about e in exactly K (their
# full stabiliser may be larger, e.g. D4 or a Klein group, without the tiling having those symmetries).
AXES = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
print("\naxial analysis (legal language):")
result = {}
for e in AXES:
    about = frozenset(R for R in frames if m.act(R, e) == e)          # rotations about e, with the identity
    half = [R for R in about if R != ident and m.act(R, m.act(R, (1, 2, 3))) == (1, 2, 3)][0]
    for K in sorted({frozenset({ident, half}), about}, key=len):
        t = group_type(list(K))
        nodes = {w for w, S in stab_all.items() if S & about == K}
        edges = {}
        for w in nodes:
            for j in (-1, 0, 1):
                o = tuple(1 + j * e[k] for k in range(3))
                if not all(0 <= x <= 2 for x in o):
                    continue
                u = window(w, o)
                if u in nodes:
                    edges.setdefault(w, []).append((u, j))
        comps = [c for c in scc(nodes, edges) if len(c) > 1 or any(u == v for v in c for u, _ in edges.get(v, ()))]
        mixed = [sorted({j for v in c for u, j in edges.get(v, ()) if u in c}) for c in comps]
        print(f"  {t:10s} axis {e}: {len(nodes)} blocks, cyclic SCCs with shift sets {mixed}")
        result.setdefault(t, False)
        result[t] |= any(len(sh) > 1 for sh in mixed)
print("axial groups occurring in the hull:", result)

# ---------------------------------------------------------------- independent check on real patches
# For each occurring group, take a mixed cycle, build the level-0 patch sigma^N(block at level N) on
# unit cells around v_0, and verify directly that every element of K maps the patch's tile poses onto
# tile poses of the patch (on the part of the patch whose image stays inside it), and that no other
# frame does.
def find_mixed_cycle(S):
    t = group_type(list(S))
    gen = [R for R in S if R != ident][0]
    e = axis_of(gen)[0]
    nodes = {w for w, T in stab_all.items() if T == S}
    edges = {}
    for w in nodes:
        for j in (-1, 0, 1):
            o = tuple(1 + j * e[k] for k in range(3))
            if all(0 <= x <= 2 for x in o) and window(w, o) in nodes:
                edges.setdefault(w, []).append((window(w, o), j))
    # BFS for a cycle using two shifts
    for start in nodes:
        stack = [(start, [], [])]
        seen = set()
        while stack:
            v, ws, js = stack.pop()
            for u, j in edges.get(v, ()):
                if u == start and len(set(js + [j])) > 1:
                    return e, ws + [v], js + [j]
                if (u, len(js)) not in seen and len(js) < 6:
                    seen.add((u, len(js)))
                    stack.append((u, ws + [v], js + [j]))
    return None


def expand(block, levels, shifts, e):
    """block at level `levels` centred on vertex 0; descend with window shifts; return level-0
    labelled cells {cell: label} of sigma^levels(block) in coordinates where v_0 = 0."""
    conf = dict(zip(POS, block))
    centre = (0, 0, 0)
    for n in range(levels):
        conf = {m.add(m.scale(2, x), d): phi[l, d] for x, l in conf.items() for d in BITS}
        centre = m.scale(2, centre)
    # v_0 position: block centre at level L is v_L; v_{n} = v_{n+1} + j_n 2^n e  (n < L)
    v = (0, 0, 0)
    for n in reversed(range(levels)):
        v = m.add(v, m.scale(shifts[n] * 2 ** n, e))
    # shift coordinates so v_0 = 0
    return {m.add(x, m.scale(-1, v)): l for x, l in conf.items()}


def poses_in(conf):
    out = set()
    for x, l in conf.items():
        f, c = labels[l]
        out.add((f, m.add(x, m.scale(-1, m.cell_lower(f, (0, 0, 0), c)))))
    return out


def rot_pose(R, pose):
    f, t = pose
    # the pose x -> f x + t maps to x -> R f x + R t ; cells are images of unit cubes, so compare cells
    return (m.compose(R, f), m.act(R, t))


print("\ndirect verification on expanded patches:")
for S in {S for S in stab_all.values() if group_type(list(S)) in ("C2(diag)", "C4")}:
    res = find_mixed_cycle(S)
    if not res:
        continue
    e, cyc, js = res
    # the cycle read backwards gives w_n = W_{j}(w_{n+1}); unroll it 3 times
    reps = max(1, 6 // len(js))
    blocks_down = []
    shifts = []
    w = cyc[0]
    seq_w = cyc * reps
    seq_j = js * reps
    # seq: seq_w[k] --j_k--> seq_w[k+1]; so the top block is seq_w[0] at level len(seq), going down
    L = len(seq_j)
    shifts_by_level = {}
    for k in range(L):
        shifts_by_level[L - 1 - k] = seq_j[k]   # level n+1 -> n with n = L-1-k
    conf = expand(seq_w[0], L, shifts_by_level, e)
    P = poses_in(conf)
    radius = 2 ** (L - 1) // 2
    inner = {p for p in P
             if all(max(abs(c) for c in m.cell_lower(p[0], p[1], x)) < radius for x in m.CHAIR_CELLS)}
    syms = [R for R in frames if all(rot_pose(R, p) in P for p in inner)]
    print(f"  {group_type(list(S)):9s} axis {e} cycle shifts {js}: tiles checked {len(inner)}, "
          f"frames preserving the patch: {group_type(syms)} (expected {group_type(list(S))})")

vertex_groups = {k[0] for k in hist}
found = {t for t, ok in result.items() if ok} | vertex_groups
expected = {"1", "C2(diag)", "C4", "D4"}
print("symmetry groups occurring:", sorted(found))
ok = found == expected and set(examples) == {"1", "D4"}
print("STATUS:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
