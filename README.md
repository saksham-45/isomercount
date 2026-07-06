# Combinatorial Enumeration: Cyclic Graphs & Snake Polyominoes

Two enumeration problems solved **exactly**, using **transfer matrices** and **Burnside's lemma** — no brute force. For each problem there is now an exact closed form, verified against 1,000 terms of ground-truth enumeration with zero mismatches.

| Problem | Object counted | Result |
|---------|----------------|--------|
| **Cyclic graphs** | Oriented cycles on `n` edges up to the dihedral group | Closed form = OEIS [A053656](https://oeis.org/A053656); computes `a(3,000,000)` (903,084 digits) in ~36 ms |
| **Linear snakes** | Edge-colored rows of `L` cells up to a 32-element symmetry group | Order-12 integer recurrence; growth rate `μ = (5+3√5)/2 ≈ 5.854` |

> **Start here:** [`CLOSED_FORM.md`](CLOSED_FORM.md) for the derivations, [`src/analysis/closed_form.py`](src/analysis/closed_form.py) for the counters, and [`docs/final_report/`](docs/final_report/) for the full written reports.

---

## Quick start

```python
from src.analysis.closed_form import cyclic_adjacency_count, linear_snake_count

cyclic_adjacency_count(8)      # 22
linear_snake_count(8)          # 653758

cyclic_adjacency_count(10**6)  # 301,024-digit integer, milliseconds
linear_snake_count(1000)       # 768-digit integer, milliseconds
```

The two core counters are **pure Python with no dependencies**. The Burnside enumerators and drawing tools additionally use NumPy and (optionally) matplotlib.

---

## Repository layout

```
src/cyclic/     oriented-cycle counting (transfer matrix, Burnside, CLI)
src/snakes/     snake-polyomino counting (Burnside, drawing, export)
src/analysis/   closed forms, recurrence discovery, symbolic verification
tests/          test suites
docs/           LaTeX sources; docs/final_report/ holds the two reports + figures
output/pdf/     compiled PDFs (reports and older proposals)
data/, figures/ generated tables and PNGs
```

Root-level scripts (`cyclic_sequences.py`, `draw_snake_graph.py`, …) are thin wrappers that re-export the implementations in `src/` for backward-compatible CLI use.

---

## The closed forms

### Cyclic graphs — a(n)

The six-state adjacency transfer matrix `M` has characteristic polynomial `χ_M(x) = x⁵(x−2)`, so `tr(Mⁿ) = 2ⁿ`. In fact the adjacency rule is **vacuous** (every configuration satisfies it), and Burnside's lemma gives

```
a(n) = (1/2n) · [ Σ_{d|n} φ(n/d)·2^d  +  (n even) · (n/2)·2^(n/2) ]
```

This is OEIS **A053656** (confirmed against the OEIS). Asymptotically `a(n) ~ 2ⁿ / (2n)`.

```
a(1..15) = 1, 2, 2, 4, 4, 9, 10, 22, 30, 62, 94, 192, 316, 623, 1096
```

### Linear snakes — N_L

For `L ≥ 2`, the sequence satisfies an order-12 integer recurrence whose characteristic polynomial factors over ℚ as

```
P(x) = (x−3)(x−1)(x²−3)(x²−5x−5)(x²−3x+1)(x⁴−5x²−5)
```

The dominant eigenvalue and amplitude are exact algebraic numbers:

```
N_L ~ A·μ^L ,   μ = (5+3√5)/2 ≈ 5.8541020 ,   A = (5+2√5)/20 ≈ 0.4736068
```

The subdominant real roots `(3±√5)/2 = φ^±2` tie the sequence to the Lucas numbers. This colored-count sequence does **not** currently appear in the OEIS (A002013 counts snake *shapes*, a different quantity).

```
N(1..15) = 4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678,
           22387074, 131052313, 767211817, 4491420695, 26293679325, 153927402355
```

### Using the counters

```python
from src.analysis.closed_form import (
    cyclic_adjacency_count,   # exact a(n), pure integer, O(√n) + one 2^n
    linear_snake_count,       # exact N_L via the order-12 recurrence
    linear_snake_sequence,    # [N_1, ..., N_max]
    linear_snake_asymptotic,  # A·μ^L to given precision (uses mpmath)
)
```

---

# Part 1: Cyclic graph sequences

Count distinct **oriented cycles** on `n` edges up to rotation and reflection. Each edge is directed (0 or 1); each vertex gets a signature in `{−2, 0, +2}` from net flow.

- **Cycle:** `n` vertices, `n` directed edges (each 0 or 1)
- **Vertex signature:** `2·(in-degree) − 2 ∈ {−2, 0, +2}`
- **Symmetry:** dihedral group `D_n` (`n` rotations + `n` reflections; reflection also flips 0↔1)

Optional constraints (all composable):

| Constraint | Description |
|------------|-------------|
| **Adjacency** | No two consecutive same-sign nonzero signatures (automatically satisfied — see the report) |
| **Primitivity** | Pattern is not a repeat of a shorter one (Möbius inversion) |
| **Forbidden subsequence** | No cyclic substring matches a canonical pattern from smaller `n` |

### Usage

```bash
# Burnside counts for n = 1..20, with optional constraints
python3 cyclic_sequences.py --max-n 20
python3 cyclic_sequences.py --max-n 20 --adjacency
python3 cyclic_sequences.py --max-n 20 --adjacency --primitive
python3 cyclic_sequences.py --max-n 20 --all-constraints

# Enumerate and draw the distinct cycles for a given n
python3 generate_cyclic_graphs.py 4 --adjacency --draw --save cyclic_n4.png
```

### Files

| File | Purpose |
|------|---------|
| `src/cyclic/cyclic_sequences.py` | Transfer matrix, Burnside count, brute enumeration, CLI |
| `src/cyclic/generate_cyclic_graphs.py` | Enumerate cycles with full labeling; draw with matplotlib |
| `src/analysis/closed_form.py` | Exact O(√n) closed-form counter |
| `tests/test_cyclic_sequences.py` | Golden-count tests; fast path vs brute |

### Performance

- Matches brute force and A053656 for small `n`.
- Closed form: `a(3,000,000)` (903,084 digits) in ~36 ms; brute force already needs >11 s at `n = 20`.
- Complexity: `O(√n)` divisor work plus one `2ⁿ` big-integer, per value.

---

# Part 2: Snake polyomino enumeration

Count distinct **edge-colored snake polyominoes**: `L` unit cells in a row (or bent chain), each cell's 4 edges (S, E, N, W) colored 0 or 1, up to a 32-element symmetry group.

- **Cell:** unit square, 4 edges, each 0 or 1; admissible if `sum ≤ 2` (11 patterns)
- **Hinge:** consecutive cells share one edge; values must match
- **Symmetry:** 8 dihedral × path reversal × global 0↔1 flip

| Mode | Shapes | L = 1..5 totals |
|------|--------|-----------------|
| **Linear** (focus of the closed form) | straight row | 4, 21, 109, 586, 3326, … |
| **Full snake** | all bent shapes (A002013 shapes) | 4, 21, 296, 3830, 54713, … |

### Usage

```bash
# Count table for L = 1..10
python3 draw_snake_graph.py --table 10 --linear
python3 draw_snake_graph.py --table 10            # full (bent) snakes

# Per-shape breakdown, drawing, JSON export
python3 draw_snake_graph.py 5 --full --linear
python3 draw_snake_graph.py 3 --draw --linear --save snakes_linear_L3.png
python3 draw_snake_graph.py 5 --linear --export linear_L5.json --limit 100
```

Common flags: `--linear`, `--full`, `--all-shapes`, `--draw`, `--save PATH`, `--export PATH`, `--limit N`, `--cols N`, `--parallel`. In drawings, **blue = 0** and **red = 1**.

### Files

| File | Purpose |
|------|---------|
| `src/snakes/draw_snake_graph.py` | Enumeration, Burnside, drawing, export |
| `src/snakes/compute_components.py` | Fast fixed-point counting for the linear arrangement |
| `src/snakes/regenerate_linear_snake_sequence.py` | Regenerate the linear sequence; verify recurrences |
| `src/analysis/closed_form.py` | Exact order-12 recurrence counter |
| `tests/test_linear_snake.py` | Brute-force vs Burnside vs transfer matrix |

### Performance

- Recurrence: `N_1000` (768 digits) in ~3 ms; `N_100000` (76,746 digits) in ~17 s — all exact integers.
- Complexity: `O(L)` big-integer steps, empirically `O(L^2.585)` with Karatsuba arithmetic.

---

# Reports

Two self-contained, publication-style write-ups live in [`docs/final_report/`](docs/final_report/) (final PDFs in `output/pdf/`):

- **`cyclic_report.tex`** — oriented cyclic graphs: transfer matrix, the `x⁵(x−2)` spectrum, Burnside reduction, the A053656 closed form, and a brute-vs-Burnside-vs-closed-form benchmark.
- **`snake_report.tex`** — linear snakes: the order-12 minimal recurrence via Berlekamp–Massey, complete root structure, exact residues, and the golden-ratio / Lucas connection.

Build them with `make cyclic_report`, `make snake_report`, or `make reports`.

---

# Tests

```bash
python3 -m tests.test_cyclic_sequences   # cyclic: golden counts, fast path vs brute
python3 -m tests.test_linear_snake       # snake: brute vs Burnside vs transfer matrix
python3 -m tests.test_closed_form        # closed forms vs Burnside (n ≤ 50, L ≤ 500)
```

---

# Dependencies

- **Python 3.7+** — the closed-form counters need nothing else.
- **NumPy** — Burnside enumerators for the snake module.
- **matplotlib** *(optional)* — `--draw` in both modules and the report figures.
- **SymPy, mpmath** *(optional)* — symbolic verification and high-precision asymptotics in `src/analysis/`.

```bash
pip install numpy matplotlib sympy mpmath
```

---

# Concepts

- **Transfer matrix** — count configurations as walks in a small state graph, via `Mⁿ`.
- **Burnside's lemma** — count orbits under symmetry by averaging fixed points.
- **Berlekamp–Massey (exact rationals)** — recover the minimal linear recurrence of a sequence without floating-point ambiguity.
- **Characteristic-polynomial factorization + Vandermonde solve** — turn a recurrence into an exact closed form with algebraic residues.

---

# References

- OEIS [A053656](https://oeis.org/A053656) — oriented cyclic graphs under the dihedral group (the cyclic sequence here)
- OEIS [A000031](https://oeis.org/A000031) — binary necklaces (used in the cyclic closed form)
- OEIS [A002013](https://oeis.org/A002013) — snake polyomino *shapes* (distinct from the colored count)
- [Burnside's lemma](https://en.wikipedia.org/wiki/Burnside%27s_lemma) · [Transfer-matrix method](https://en.wikipedia.org/wiki/Transfer-matrix_method)
- Pisanski, Schattschneider, Servatius, *Applying Burnside's lemma to a one-dimensional Escher problem*, Math. Mag. **79** (2006), 167–180.
