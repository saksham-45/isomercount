# Isomer count — oriented cycles under dihedral symmetry

I’m counting distinct oriented cycles on **n** edges up to dihedral symmetry (rotation + reflection). You can turn on extra constraints: **adjacency** (no two consecutive same-sign nonzero vertices), **primitivity**, and **forbidden subsequences**. The main sequence I care about is the adjacency-only one, and I compute it in polynomial time using a transfer matrix and Burnside’s lemma—no brute-force enumeration for large n.

---

## What I got

- **Small n:** I checked a(1..15) against brute force and known values (1, 2, 2, 4, 4, 9, 10, 22, 30, 62, …). The optimized path matches.
- **Big n:**
  - **a(1,000,000)** — ran in a couple of seconds (hundreds of thousands of digits).
  - **a(2,000,000)** — **602,054 digits**, about **2.8 seconds**.
  - **a(3,000,000)** — **903,084 digits**, about **4.5 seconds**.
- I didn’t use any closed formula for the main count—everything comes from the combinatorial setup. Complexity for a single n is roughly **O(d(n)·√n)** time and **O(d(n))** space (divisors + totients + a tiny matrix).

---

## What I’m actually counting

Think of a cycle with **n** vertices and **n** directed edges. Each edge is either “forward” (0) or “backward” (1). At each vertex I get a number in **{−2, 0, +2}** from the two edges meeting there (net in vs out). Valid cycles are automatically balanced. I count these things up to **rotation** and **reflection** (reverse the cycle and flip every edge).

The constraints I added:
1. **Adjacency** — no two consecutive nonzero vertex values with the same sign.
2. **Primitivity** — the cycle isn’t just a shorter pattern repeated (I use Möbius inversion for that).
3. **Forbidden subsequence** — no cyclic substring matches a canonical pattern from a smaller length (for recursive constraints).

The code is built around the **adjacency** case and the fast (matrix + Burnside) path.

---

## Why I did it this way

### Why not brute force?

Brute force would mean looping over all **2^n** edge choices, then reducing by symmetry. That’s **O(2^n · n)**. Fine for n ≈ 15, hopeless for millions. So I wanted a way to **count** without listing everything.

### Why a transfer matrix?

Along the cycle, I only need to remember a small amount: the last vertex value and the current edge. That’s **6 states** (3 vertex values × 2 edges). The adjacency rule is local—it only tells you which next states are allowed. So I built a **6×6** matrix **M** where M[s][t] = 1 if you can go from state s to state t in one step and still satisfy adjacency. Then the number of closed walks of length n is related to **trace(M^n)** (with the right closure condition). Computing **M^n** takes **O(log n)** matrix multiplies, so that part is cheap.

### Why Burnside?

I don’t want “closed sequences”—I want **orbits** under the dihedral group (n rotations + n reflections). Burnside says: number of orbits = (1/|G|) × (sum over group elements of how many configurations they fix). For **rotations**: a config is fixed by rotation by k iff it has period **gcd(n,k)**. So I need the count of closed sequences of length **d** for each divisor **d** of **n**, weighted by **φ(n/d)**. I get those counts from the matrix. For **reflections** I count sequences fixed by “reverse + flip” using one more matrix power (length n/2) and a boundary check. So the final count is **(rot_sum + reflection_term) / (2n)**. No enumeration—just a loop over divisors, a few matrix powers, and some totients.

### Why O(√n) for divisors and totient?

At first I had **divisors(n)** as “loop 1 to n and check n % d == 0”—that’s **O(n)**. For n = 2 million that’s 2 million steps. So I changed it to loop only up to **√n** and for each divisor d I also get n/d. Same list, **O(√n)** time.

For **φ(m)** I was doing “count k in 1..m with gcd(k,m)=1,” which is **O(m)**. I switched to the formula **φ(m) = m · Π(1 − 1/p)** over primes p dividing m, and I get those primes by trial division up to **√m**. So **O(√m)** per totient. With that, a single run for n = 2e6 dropped to a few seconds instead of choking on millions of steps.

### Why no formula for the main count?

A053656 (no constraints) has a known formula. I deliberately don’t use it for the main sequence—I wanted the program to derive the count purely from the combinatorial construction (edges → vertex signatures → symmetry). I kept the formula in the repo only for **verification** (e.g. checking a(1..15) when there are no constraints).

---

## What’s in the repo

- **`cyclic_sequences.py`** — Main code: transfer matrix, Burnside count, brute enumeration for small n, constraints, and the CLI.
- **`a053656_binary.py`** — Reference for A053656: binary encoding, dihedral canonical form, and the formula (verification only).
- **`test_cyclic_sequences.py`** — Tests with golden counts for small n and different constraint combos; checks that the fast path matches brute where we can.
- **`run_custom.py`** — Small script to run the counter (e.g. for one large n).

I didn’t commit the huge output files (full decimals for a(2e6), a(3e6)); the digit counts and runtimes are in this README.

---

## How to run

```bash
# Counts for n=1..max_n with optional constraints
python3 cyclic_sequences.py --max-n 20 [--adjacency] [--primitive] [--forbidden] [--verbose]

# For one huge n from Python:
# from cyclic_sequences import count_adjacency_burnside
# count_adjacency_burnside(2_000_000)
```

Tests (no pytest needed):

```bash
python3 test_cyclic_sequences.py
```

---

## Complexity in a sentence

For a single **n**: time **O(d(n)·√n + d(n)·log n)** (divisors, totients, matrix), space **O(d(n))**; the answer a(n) itself has about **Θ(n)** digits if you write it out.

That’s the project—why it’s fast and why I built it this way.
