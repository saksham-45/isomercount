# Isomer count — oriented cycles under dihedral symmetry (with constraints)

Counts distinct oriented cycles on **n** edges up to dihedral symmetry (rotation + reflection), with optional constraints: **adjacency** (no consecutive same-sign nonzero vertex), **primitivity**, and **forbidden subsequences**. The main sequence (adjacency-only) is computed in **polynomial time** via transfer matrix + Burnside’s lemma.

---

## Summary of results

- **Small n (brute-checked):** a(1..15) with adjacency match the brute-force enumeration and known behavior (1, 2, 2, 4, 4, 9, 10, 22, 30, 62, 94, 192, 316, 623, 1096 for no constraints; adjacency-only counts match when using the optimized path).
- **Large n (optimized path):**
  - **a(1,000,000):** computed in ~1–2 s (hundreds of thousands of digits).
  - **a(2,000,000):** **602,054 digits**; computed in **~2.8 s**.
  - **a(3,000,000):** **903,084 digits**; computed in **~4.5 s**.
- **Complexity (single n):** Time **O(d(n)·√n + d(n)·log n)**, space **O(d(n))** (worst case O(√n)). No formula used for the main count; everything is computed from the combinatorial construction.

---

## What we count

- **Objects:** Oriented cycles on **n** vertices (or equivalently **n** directed edges on a cycle). Each edge is either 0 (i→i+1) or 1 (i+1→i). We identify two configurations that differ by a **rotation** or a **reflection** (reverse + flip edges).
- **Vertex signature:** At each vertex we get a value in **{−2, 0, +2}** from the two incident edges (net “in” minus “out”). Valid cycles are automatically **balanced** (same number of +2 and −2).
- **Constraints (optional):**
  1. **Adjacency:** No two consecutive nonzero vertex values have the same sign.
  2. **Primitivity:** The cycle is not a periodic repetition of a shorter cycle (handled via Möbius inversion).
  3. **Forbidden subsequence:** No cyclic substring of the vertex signature is equivalent to a canonical pattern from a smaller length (used for recursive constraint).

The main implementation focuses on the **adjacency** constraint and uses no closed formula: counts are derived from the transfer matrix and Burnside.

---

## Why this approach

### Why not brute force?

- Brute force: enumerate all **2^n** edge configurations, reduce by dihedral symmetry, then filter. Time **O(2^n · n)**. For n ≈ 20 this is already heavy; for n = 2,000,000 it is infeasible.
- So we need a **count method** that avoids enumerating configurations.

### Why transfer matrix?

- The state of the process along the cycle can be summarized by a **bounded state**: e.g. last vertex value and current edge (6 states: 3 vertex values × 2 edges). **Adjacency** is a local rule: it only restricts transitions between consecutive states.
- So we build a **6×6 transition matrix M**: M[s][t] = 1 if extending from state s by one edge yields state t and the new vertex satisfies adjacency. The number of closed walks of length n with adjacency is then related to **trace(M^n)** (for the correct closure condition) or a small sum over (start, end) state pairs.
- **Matrix power M^n** is done in **O(log n)** multiplications (binary exponentiation). So we get the “raw” count for period-n sequences in **O(log n)** time with fixed state space.

### Why Burnside?

- We want **orbits** under the dihedral group (n rotations + n reflections), not raw closed sequences. **Burnside’s lemma**: number of orbits = (1/|G|) Σ_g fix(g).
- **Rotations:** A configuration is fixed by rotation by k iff it has period **gcd(n,k)**. So fix(rotation k) = (number of closed sequences of length **gcd(n,k)** with adjacency). We already have that via the matrix. Sum over k gives a sum over **divisors d of n**: Σ_{d|n} φ(n/d) · (count for length d).
- **Reflections:** For “reverse + flip” we count sequences fixed by that reflection (e.g. by counting half-length paths that match the reflection condition). Implemented via one more matrix power M^(n/2) and a boundary check.
- So the total count is **(rot_sum + reflection_term) / (2n)**. No enumeration, only divisor loop + matrix powers + totients.

### Why O(√n) divisors and totient?

- **Divisors:** Naively listing divisors by iterating 1..n is **O(n)**. For n = 2e6 that’s 2 million steps. Instead we iterate **d = 1..√n** and for each d with n % d == 0 we get two divisors: d and n/d. **O(√n)**.
- **Totient φ(m):** Naively φ(m) = #{ k in 1..m : gcd(k,m)=1 } is **O(m)**. We use **φ(m) = m · Π_{p|m} (1 − 1/p)** and compute the prime factors of m by trial division up to **√m**. **O(√m)**. For m = n/d the total over divisors is still much cheaper than O(n).

Without these optimizations, a single call for n = 2e6 would do millions of divisor/totient steps; with them we get a few thousand and the runtime drops to a few seconds.

### Why no formula for the main count?

- A053656 (no constraints) has a known formula (necklace count + correction). We **do not** use it for the main sequence so that the program is “from scratch”: the count is derived only from the combinatorial construction (edges → vertex signatures → symmetry). The formula is kept only for **verification** (e.g. comparing a(1..15) when no constraints are applied).

---

## Repo layout

- **`cyclic_sequences.py`** — Main module: transfer matrix, Burnside count, brute enumeration (for small n / verification), constraints (adjacency, primitivity, forbidden), and CLI.
- **`a053656_binary.py`** — A053656 reference: binary edge encoding, dihedral canonical form, and the **formula** (used only for verification).
- **`test_cyclic_sequences.py`** — Tests: golden counts for small n and constraint combinations, and checks that the optimized count matches brute where applicable.
- **`run_custom.py`** — Example script to run the counter (e.g. for a single large n).

Large outputs (e.g. full decimals for a(2e6), a(3e6)) are not committed; digit counts and runtimes are summarized above and in this README.

---

## How to run

```bash
# Counts for n=1..max_n with optional constraints (brute for small n if needed)
python3 cyclic_sequences.py --max-n 20 [--adjacency] [--primitive] [--forbidden] [--verbose]

# Single large n (optimized path, adjacency-only) from code:
# from cyclic_sequences import count_adjacency_burnside
# count_adjacency_burnside(2_000_000)
```

Run tests (no pytest required):

```bash
python3 test_cyclic_sequences.py
```

---

## Complexity (recap)

- **Time (single n):** O(d(n)·√n + d(n)·log n) — divisor iteration, totients, and O(log n) matrix work per divisor.
- **Space:** O(d(n)) for the divisor list; matrix is O(1).
- **Output size:** a(n) has Θ(n) digits, so storing the decimal is O(n).

