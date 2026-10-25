# Hull and symmetry of Chair44 tilings

These scripts give computer-assisted answers to Open Problems 1 and 2 of arXiv:2609.19214 (§10):

- **Hull.** Every tiling by Q lies in the substitution hull. So every Chair44 tiling is pure point diffractive,
  and the space of tilings is minimal and uniquely ergodic.
- **Symmetry.** The symmetry group of a tiling by Q is conjugate to exactly one of 1, C2 (about a face-diagonal axis), C4
  and D4, and each occurs. So |Sym(T)| ≤ 8.

The scripts use only the Python standard library and this repository's own data:
- `certificates/candidate_certificate.json` (the 44-contact atlas and the children);
- `simulations/explore/framed_chair.py`;
- `verify/substitution_modular_coincidence.py`.

```
python3 verify/hull/chair44_vertex_blocks.py   blocks.json   # admissible 2x2x2 vertex blocks: 360189   (~8 min)
python3 verify/hull/chair44_hull_prune.py      blocks.json   # Win(sigma(admissible)) = the 1410 legal blocks
python3 verify/hull/chair44_symmetry_census.py               # groups 1, C2(diag), C4, D4; checked on real patches
```

Each script ends with `STATUS: PASS` or `FAIL`.

## Argument

1. A *block* is the labelling of the 8 unit cells at a lattice vertex.
   - It is *admissible* if the chairs it names are consistent and disjoint, and all face contacts are in A44.
     By Thm 5.6 and Lemma 5.3, every block of a tiling by Q is admissible.
   - It is *legal* if it occurs in some σⁿ(j).
2. **Lemma (computed).** Every window of σ(b), for b admissible, is legal.
3. **Hull.** By Thms 6.4–6.5, T = σ(D(T)) with D(T) again a tiling by Q.
   - Every patch of T lies in σⁿ(β) for a block β of Dⁿ(T).
   - β is a window of σ(β′) with β′ admissible, so β is legal by the lemma.
   - Hence every patch of T is legal.
4. **Symmetry.**
   - No element of order 3: its axis meets a tile interior, and Q has no self-isometry. So Sym(T) ≤ D4.
   - Point groups fix a vertex of every supertile lattice. The blocks there are periodic under the central window, and those blocks are the
     27 fixed seeds, with stabilisers 1 and D4.
   - Cyclic groups have an axis, with direction e, through vertices of every level. The group K occurs iff a certain graph has a
     strongly connected component using two different axial window shifts. The graph's vertices are the legal blocks whose
     stabiliser meets the rotations about e in exactly K; their full stabiliser may be larger. That holds for C2(diag) and C4,
     but not for C2 about a coordinate axis.
   - Independent check: Brent Hartshorn's `spectre` repository regenerates the tile and its 44-contact atlas from scratch
     (`symchair44.py`). From that atlas, `chair44_hull_admissible.py` gets the same 360,189 admissible blocks, and all windows
     of their substitutions are legal.

Paper draft: T. Romanov, *Every Chair44 tiling is a substitution tiling, and its symmetry groups* (2026). It is
computer-assisted and unrefereed, and independent checks are welcome.
