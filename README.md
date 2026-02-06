# Isomer Count: Counting Oriented Cycles with Dihedral Symmetry

Hi! This project is a specialized tool I built to count distinct oriented cycles on $n$ edges. Specifically, it handles the tricky math of **dihedral symmetry** (accounting for both rotations and reflections) while applying various constraints like adjacency rules and primitivity.

The core of the project is an optimized "path" that uses a Transfer Matrix approach combined with Burnside's Lemma. This allows me to jump from simple brute-force checks for small $n$ to calculating results for $n = 3,000,000$ in just a few seconds.

---

## 🚀 Performance Highlights

I’ve optimized the math so that even for massive values of $n$, the computation stays fast:

* **Small $n$:** Verified against brute-force enumeration for $a(1..15)$.
* **Large $n$:** * **$n = 1,000,000$:** Computed in ~1–2 seconds.
    * **$n = 3,000,000$:** Generated a result with **903,084 digits** in about 4.5 seconds.
* **Complexity:** The time complexity is $O(d(n)\sqrt{n} + d(n)\log n)$, where $d(n)$ is the number of divisors. It’s built to be lean and efficient.

---

## 🔍 What am I counting?

I’m looking at oriented cycles with $n$ vertices (or $n$ directed edges). Each edge is either 0 or 1. I want to find the number of unique configurations after we ignore versions that are just rotations or reflections (flipping and reversing) of each other.

### The Constraints
I’ve added support for three specific types of constraints:
1.  **Adjacency:** No two consecutive non-zero vertex values can have the same sign.
2.  **Primitivity:** I filter out cycles that are just periodic repetitions of shorter cycles (using Möbius inversion).
3.  **Forbidden Subsequences:** I can block specific cyclic patterns from appearing in the vertex signature.

---

## 🛠️ My Approach

### Why not just brute force it?
Brute force is fine for $n=20$, but at $n=2,000,000$, it’s impossible. Enumerating $2^n$ configurations would take lifetimes. I needed a way to count without actually listing every possibility.

### The Transfer Matrix Trick
Instead of listing cycles, I use a **6x6 transition matrix ($M$)**. By looking at the states of the vertices and edges, I can turn the counting problem into a matrix power problem. Using binary exponentiation, I can find the "raw" count in $O(\log n)$ time.

### Burnside’s Lemma
To handle the symmetry, I use **Burnside's Lemma**. 
* **For Rotations:** I sum the counts for lengths based on the divisors of $n$.
* **For Reflections:** I use a separate matrix power $M^{n/2}$ to count configurations that stay the same when flipped.

By combining these, I get the total orbit count using the formula:
$$\frac{\text{rotation\_sum} + \text{reflection\_term}}{2n}$$

---

## Getting Started
The main implementation focuses on the **adjacency** constraint. Since I'm using a combinatorial construction rather than a closed-form formula, the code is flexible enough to handle different constraints without a total rewrite.
