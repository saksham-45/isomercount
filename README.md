# Combinatorial Enumeration: Cyclic Graphs & Snake Polyominoes

Two complementary enumeration projects using **transfer matrices** and **Burnside's lemma** to count distinct structures under symmetry—no brute force, verified for correctness.

## Repository layout

- `src/cyclic/`: cyclic (ring) counting implementation
- `src/snakes/`: snake/linear counting implementation
- `src/analysis/`: recurrence/asymptotic/verification helpers
- `tests/`: test modules
- `docs/`: LaTeX sources and bibliography
- `figures/`: generated PNGs
- `data/`: generated tables/exports/large outputs
- `output/pdf/`: compiled PDFs (`ass1.pdf`, `ass2.pdf`, etc.)

Root-level scripts like `cyclic_sequences.py` and `draw_snake_graph.py` are thin wrappers kept for backwards-compatible CLI usage.

---

## Overview

| Project | What it counts | Key features |
|---------|----------------|--------------|
| **Cyclic graphs** | Oriented cycles on n edges under dihedral symmetry | Transfer matrix, O(√n) divisors, handles n = millions |
| **Snake graphs** | Edge-colored chains of squares (polyominoes) | Linear & bent shapes, drawing, export, Burnside orbits |

---

# Part 1: Cyclic Graph Sequences

Count distinct **oriented cycles** on n edges up to rotation and reflection. Each edge is directed (0 or 1); vertices get signatures in {−2, 0, +2} from net flow.

## What we count

- **Cycle**: n vertices, n directed edges (each 0 or 1)
- **Vertex signature**: at each vertex, count incoming edges → −2, 0, or +2
- **Symmetry**: dihedral group (n rotations + n reflections)
- **Valid**: automatically balanced; optional constraints below

## Constraints

| Constraint | Description |
|------------|-------------|
| **Adjacency** | No two consecutive nonzero vertex values with the same sign |
| **Primitivity** | Pattern is not a repeat of a shorter pattern (Möbius inversion) |
| **Forbidden subsequence** | No cyclic substring matches a canonical pattern from smaller n |

## Usage

```bash
# Counts for n=1..20 with optional constraints
python3 cyclic_sequences.py --max-n 20
python3 cyclic_sequences.py --max-n 20 --adjacency
python3 cyclic_sequences.py --max-n 20 --adjacency --primitive
python3 cyclic_sequences.py --max-n 20 --all-constraints
python3 cyclic_sequences.py --max-n 20 --adjacency --verbose

# Huge n (from Python)
# from cyclic_sequences import count_adjacency_burnside
# count_adjacency_burnside(2_000_000)  # ~602k digits, ~3 seconds
```

## Files

| File | Purpose |
|------|---------|
| `src/cyclic/cyclic_sequences.py` | Main: transfer matrix, Burnside count, brute enumeration, CLI |
| `src/cyclic/generate_cyclic_graphs.py` | Enumerate cycles with full labeling; draw with matplotlib |
| `src/cyclic/a053656_binary.py` | A053656 formula reference (verification only) |
| `tests/test_cyclic_sequences.py` | Golden-count tests; fast path vs brute |

## Performance

- **Small n:** Matches brute force and A053656 for n ≤ 15
- **Large n:** a(2,000,000) ≈ 602k digits in ~3 seconds
- **Complexity:** O(d(n)·√n + d(n)·log n) per n; space O(d(n))

## Generate & draw cyclic graphs

```bash
# List all unique cycles for n=4
python3 generate_cyclic_graphs.py 4

# With adjacency constraint
python3 generate_cyclic_graphs.py 6 --adjacency

# With adjacency + primitivity
python3 generate_cyclic_graphs.py 5 --adjacency --primitive

# Draw with matplotlib
python3 generate_cyclic_graphs.py 4 --adjacency --draw --save cyclic_graphs_n4.png

# Write to file
python3 generate_cyclic_graphs.py 5 --adjacency -o graphs_n5.txt
```

---

# Part 2: Snake Graph Enumeration

Count distinct **edge-colored snake polyominoes**: L square cells arranged in chains or rows, each cell with 4 edges (S,E,N,W) colored 0 or 1. Equivalence under 32-element symmetry (dihedral × reversal × flip).

## What we count

- **Cell**: unit square with 4 edges (S,E,N,W), each 0 or 1
- **Valid patterns**: sum(edges) ≤ 2 per cell (11 patterns)
- **Hinge**: consecutive cells share exactly one edge; values must match
- **Symmetry**: 8 dihedral × path reversal × global 0↔1 flip

## Modes

| Mode | Shapes | L=1..10 totals |
|------|--------|----------------|
| **Full snake polyomino** | All bent shapes (A002013) | 4, 21, 296, 3,830, 54,713, … |
| **Linear** | Only straight line (L cells in a row) | 4, 21, 109, 586, 3,326, … |

## Usage

```bash
# Table: count for L=1..10
python3 draw_snake_graph.py --table 10
python3 draw_snake_graph.py --table 10 --linear

# Count shapes for a given L
python3 draw_snake_graph.py 5
python3 draw_snake_graph.py 5 --linear

# Full per-shape breakdown
python3 draw_snake_graph.py 5 --full
python3 draw_snake_graph.py 5 --full --linear

# Draw all unique snakes
python3 draw_snake_graph.py 3 --draw --save snakes_L3.png
python3 draw_snake_graph.py 3 --draw --linear --save snakes_linear_L3.png

# Export examples to JSON
python3 draw_snake_graph.py 5 --linear --export linear_L5.json --limit 100

# Draw first 100 (when total is large)
python3 draw_snake_graph.py 5 --draw --linear --limit 100 --save examples.png
```

## Options

| Flag | Description |
|------|-------------|
| `--table N` | Print count table for L=1..N |
| `--linear` | Linear only: 1 shape per L |
| `--full` | Per-shape and total distinct snakes |
| `--all-shapes` | Include non-snake shapes (U-bends etc.) |
| `--draw` | Draw all unique snakes with matplotlib |
| `--save PATH` | Save drawing to file |
| `--export PATH` | Export examples to JSON |
| `--limit N` | Draw/export only first N snakes |
| `--cols N` | Columns in drawing grid (default 5) |
| `--parallel` | Parallel processing for --table |

## Files

| File | Purpose |
|------|---------|
| `src/snakes/draw_snake_graph.py` | Main: enumeration, Burnside, drawing, export |
| `src/snakes/compute_components.py` | Fast fixed-point counting for linear arrangement |
| `src/snakes/regenerate_linear_snake_sequence.py` | Regenerate linear sequence; verify recurrences |
| `tests/test_linear_snake.py` | Verification: brute-force vs Burnside, transfer matrix |

## Colors in drawings

- **Blue** = edge value 0  
- **Red**  = edge value 1  

## Linear sequence (L=1..15)

```
4, 21, 109, 586, 3326, 19209, 111871, 653758, 3824678, 22387074, 131052313, 767211817, 4491420695, 26293679325, 153927402355
```

---

---

# Tests

```bash
python3 test_cyclic_sequences.py      # wrapper (runs tests/test_cyclic_sequences.py)
python3 test_linear_snake.py          # wrapper (runs tests/test_linear_snake.py)
```

---

# Dependencies

- Python 3.7+
- **matplotlib** (optional, for `--draw` in both modules)

```bash
pip install matplotlib
# or use the project venv: source .venv/bin/activate
```

---

# Concepts

Both projects use:

- **Transfer matrix** — local state transitions; count without listing
- **Burnside's lemma** — orbits under symmetry; avoid overcounting
- **Canonical form** — one representative per equivalence class
- **Constraints** — local rules (adjacency, forbidden patterns) that compose

---

# References

- OEIS [A053656](https://oeis.org/A053656) — oriented cycles (no constraints)
- OEIS [A002013](https://oeis.org/A002013) — snake polyomino shapes
- [Burnside's lemma](https://en.wikipedia.org/wiki/Burnside%27s_lemma)
- [Transfer matrix method](https://en.wikipedia.org/wiki/Transfer-matrix_method)
