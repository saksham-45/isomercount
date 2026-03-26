#!/usr/bin/env python3
"""
Generate all unique cyclic graphs for a given n, with full labeling.

For research-accurate enumeration of oriented cycles on n edges under dihedral
symmetry (rotation + reflection), with optional constraints: adjacency,
primitivity, forbidden subsequences.

Each graph is presented with:
  - Canonical edge sequence (0/1 encoding)
  - Per-edge direction: vertex tail → vertex head (out from tail, in to head)
  - Vertex signatures (-2, 0, +2) indicating net flow at each vertex
  - Optional: adjacency/primitivity/forbidden validation

Edge convention (from cyclic_sequences.py):
  - Edge i connects vertex i and vertex (i+1) mod n.
  - edge[i]=0: direction i → (i+1), i.e. OUT from vertex i, IN to vertex (i+1)
  - edge[i]=1: direction (i+1) → i, i.e. IN to vertex i, OUT from vertex (i+1)

Usage:
    python generate_cyclic_graphs.py 4
    python generate_cyclic_graphs.py 6 --adjacency
    python generate_cyclic_graphs.py 5 --adjacency --primitive --output graphs_n5.txt
"""

import argparse
import math
import sys
from typing import List, Tuple, Sequence, Set, Dict, Optional

# Import from the cyclic_sequences module (package-relative)
from .cyclic_sequences import (
    generate_all_oriented_cycles,
    canonical_form_sig,
    is_adjacency_valid,
    is_primitive,
    contains_forbidden,
    count_adjacency_burnside,
    count_adjacency_primitive_burnside,
)


def get_edge_direction(edge_val: int, edge_idx: int, n: int) -> Tuple[int, int, str, str]:
    """
    Return (tail_vertex, head_vertex, tail_label, head_label) for edge_idx.
    Edge i is between vertex i and vertex (i+1) mod n.
    - edge=0: i → (i+1), so tail=i, head=(i+1)
    - edge=1: (i+1) → i, so tail=(i+1), head=i
    """
    v_lo = edge_idx
    v_hi = (edge_idx + 1) % n
    if edge_val == 0:
        # Direction: v_lo → v_hi
        return (v_lo, v_hi, "out", "in")
    else:
        # Direction: v_hi → v_lo
        return (v_hi, v_lo, "out", "in")


def format_vertex_signature(sig: Sequence[int]) -> str:
    """Format vertex signature as (+2, 0, -2, ...)."""
    parts = []
    for s in sig:
        if s == 0:
            parts.append("0")
        elif s > 0:
            parts.append(f"+{s}")
        else:
            parts.append(str(s))
    return "(" + ", ".join(parts) + ")"


def format_edge_sequence(edge_seq: Sequence[int]) -> str:
    """Format edge sequence as binary string and as tuple."""
    binary = "".join(str(b) for b in edge_seq)
    tuple_str = "(" + ", ".join(str(b) for b in edge_seq) + ")"
    return binary, tuple_str


def describe_graph(
    edge_seq: Tuple[int, ...],
    vertex_sig: Tuple[int, ...],
    graph_num: int,
    total: int,
    compact: bool = False,
) -> str:
    """Produce a human-readable description of one cyclic graph."""
    n = len(edge_seq)
    lines: List[str] = []

    binary, tuple_str = format_edge_sequence(edge_seq)
    sig_str = format_vertex_signature(vertex_sig)

    header = f"Graph {graph_num}/{total}  (n={n})"
    lines.append("=" * (len(header) + 4))
    lines.append(f"  {header}")
    lines.append("=" * (len(header) + 4))

    lines.append("")
    lines.append("  Edge sequence (0=→, 1=←):")
    lines.append(f"    Binary: {binary}")
    lines.append(f"    Tuple:  {tuple_str}")
    lines.append("")

    if not compact:
        lines.append("  Edge directions (tail → head; 'out' from tail, 'in' to head):")
        for i in range(n):
            ev = edge_seq[i]
            tail, head, tail_lbl, head_lbl = get_edge_direction(ev, i, n)
            # Always show tail → head (direction of the directed edge)
            lines.append(
                f"    Edge {i:2d}: vertex {tail} → vertex {head}  "
                f"({tail_lbl} from v{tail}, {head_lbl} to v{head})"
            )
        lines.append("")

    lines.append("  Vertex signatures (net flow: -2=all out, 0=balanced, +2=all in):")
    lines.append(f"    {sig_str}")
    for v in range(n):
        s = vertex_sig[v]
        if s == -2:
            desc = "all out"
        elif s == 0:
            desc = "balanced"
        else:
            desc = "all in"
        lines.append(f"    V{v}: {s:+d} ({desc})")
    lines.append("")

    return "\n".join(lines)


def draw_cyclic_graphs(
    cycles_list: List[Tuple[Tuple[int, ...], Tuple[int, ...]]],
    n: int,
    save_path: Optional[str] = None,
) -> None:
    """Draw all cyclic graphs using matplotlib (vertices on circle, directed edges as arrows)."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrowPatch
    except ImportError:
        print("matplotlib required for --draw. Install with: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    n_graphs = len(cycles_list)
    ncols = min(n_graphs, 4) if n_graphs <= 4 else 3
    nrows = (n_graphs + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    # Vertex positions on unit circle (vertex 0 at top, going CW)
    radius = 1.0
    positions = []
    for i in range(n):
        angle = -2 * math.pi * i / n + math.pi / 2  # 0 at top
        positions.append((radius * math.cos(angle), radius * math.sin(angle)))

    for idx, (edge_seq, vertex_sig) in enumerate(cycles_list):
        ax = axes_flat[idx]
        ax.set_aspect("equal")
        ax.axis("off")

        binary = "".join(str(b) for b in edge_seq)
        sig_str = format_vertex_signature(vertex_sig)
        ax.set_title(f"Graph {idx + 1}: {binary}\n{sig_str}", fontsize=10)

        # Draw vertices
        for v, (x, y) in enumerate(positions):
            color = "#2ecc71" if vertex_sig[v] == 0 else "#3498db" if vertex_sig[v] > 0 else "#e74c3c"
            ax.plot(x, y, "o", markersize=14, color=color, zorder=3)
            ax.text(x, y, str(v), ha="center", va="center", fontsize=11, fontweight="bold", zorder=4)

        # Draw directed edges (use curved arrows for edges that go "around" the circle)
        for i in range(n):
            tail, head, _, _ = get_edge_direction(edge_seq[i], i, n)
            x1, y1 = positions[tail]
            x2, y2 = positions[head]

            # Curvature for chord (avoid overlap when tail and head are far)
            dx, dy = x2 - x1, y2 - y1
            dist = math.sqrt(dx * dx + dy * dy)
            curve = 0.15 if dist < 1.5 else 0

            arrow = FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                arrowstyle="->,head_width=0.25,head_length=0.4",
                color="#34495e",
                linewidth=2,
                connectionstyle=f"arc3,rad={curve}",
                zorder=1,
            )
            ax.add_patch(arrow)

        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)

    for j in range(n_graphs, len(axes_flat)):
        axes_flat[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=120, bbox_inches="tight")
        print(f"Saved to {save_path}")
    else:
        plt.show()


def collect_unique_cycles(
    n: int,
    use_adjacency: bool = False,
    use_primitivity: bool = False,
    use_forbidden: bool = False,
    all_sets_so_far: Optional[Dict[int, Set[Sequence]]] = None,
) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """
    Return a list of (edge_seq, vertex_sig) representatives, one per distinct
    cyclic graph under the chosen constraints. Each is in canonical form.
    """
    if all_sets_so_far is None:
        all_sets_so_far = {}

    cycles = generate_all_oriented_cycles(n)
    seen_sigs: Set[Tuple[int, ...]] = set()
    result: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []

    for edge_seq, vertex_sig in sorted(cycles):
        canon_sig = canonical_form_sig(vertex_sig)
        if canon_sig in seen_sigs:
            continue
        if use_adjacency and not is_adjacency_valid(canon_sig):
            continue
        if use_primitivity and not is_primitive(canon_sig):
            continue
        if use_forbidden and contains_forbidden(canon_sig, all_sets_so_far):
            continue
        seen_sigs.add(canon_sig)
        result.append((edge_seq, vertex_sig))

    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate all unique cyclic graphs for given n with proper edge/vertex labeling.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Edge convention:
  - Each edge connects two consecutive vertices on the cycle (with wrap-around).
  - 0: direction from lower index to higher (vertex i → vertex i+1)
  - 1: direction from higher index to lower (vertex i+1 → vertex i)
  - "out" = edge emanates from vertex; "in" = edge enters vertex.

Vertex signature:
  - At each vertex, count how many of the two incident edges point IN.
  - 0 in → -2, 1 in → 0, 2 in → +2.
        """,
    )
    parser.add_argument(
        "n",
        type=int,
        help="Number of edges (and vertices) in the cycle",
    )
    parser.add_argument(
        "--adjacency",
        action="store_true",
        help="Apply adjacency constraint (no two consecutive same-sign nonzero vertices)",
    )
    parser.add_argument(
        "--primitive",
        action="store_true",
        help="Apply primitivity constraint (no repeating block)",
    )
    parser.add_argument(
        "--forbidden",
        action="store_true",
        help="Apply forbidden subsequence constraint",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Write output to file instead of stdout",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Compact output (skip per-edge direction listing)",
    )
    parser.add_argument(
        "--max-n",
        type=int,
        default=14,
        help="Maximum n for which full enumeration is attempted (default 14)",
    )
    parser.add_argument(
        "--draw",
        action="store_true",
        help="Draw graphs with matplotlib (vertices on circle, arrows for edge direction)",
    )
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        metavar="PATH",
        help="Save drawing to file (use with --draw)",
    )
    args = parser.parse_args()

    n = args.n
    if n < 1:
        print("Error: n must be >= 1", file=sys.stderr)
        return 1
    if n > args.max_n:
        print(
            f"Warning: n={n} exceeds --max-n={args.max_n}. "
            "Enumeration may be slow (O(2^n)). Use --max-n to override.",
            file=sys.stderr,
        )

    # Build all_sets_so_far for forbidden constraint (valid canonical sigs for 1..n-1)
    all_sets: Dict[int, Set[Tuple[int, ...]]] = {}
    if args.forbidden:
        for k in range(1, n):
            all_sets[k] = set()
            for edge_seq, vertex_sig in generate_all_oriented_cycles(k):
                canon_sig = canonical_form_sig(vertex_sig)
                if args.adjacency and not is_adjacency_valid(canon_sig):
                    continue
                if args.primitive and not is_primitive(canon_sig):
                    continue
                if k >= 2 and contains_forbidden(canon_sig, all_sets):
                    continue
                all_sets[k].add(canon_sig)

    cycles_list = collect_unique_cycles(
        n,
        use_adjacency=args.adjacency,
        use_primitivity=args.primitive,
        use_forbidden=args.forbidden,
        all_sets_so_far=all_sets,
    )

    # Expected count from fast path (when applicable)
    if args.adjacency and not args.forbidden:
        if args.primitive:
            expected = count_adjacency_primitive_burnside(n)
        else:
            expected = count_adjacency_burnside(n)
        if len(cycles_list) != expected:
            print(
                f"Warning: enumerated {len(cycles_list)} but Burnside count is {expected}",
                file=sys.stderr,
            )

    if args.draw:
        draw_cyclic_graphs(cycles_list, n, save_path=args.save)

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write(
            f"Unique cyclic graphs for n={n}\n"
            f"Constraints: adjacency={args.adjacency}, primitive={args.primitive}, forbidden={args.forbidden}\n"
            f"Total count: {len(cycles_list)}\n\n"
        )
        # Summary table: graph# | binary | vertex signature
        out.write("Summary (quick reference):\n")
        out.write(f"  {'#':>3}  {'Edge (binary)':>12}  Vertex signature\n")
        out.write("  " + "-" * 50 + "\n")
        for i, (edge_seq, vertex_sig) in enumerate(cycles_list, 1):
            binary = "".join(str(b) for b in edge_seq)
            sig_str = format_vertex_signature(vertex_sig)
            out.write(f"  {i:>3}  {binary:>12}  {sig_str}\n")
        out.write("\n")
        for i, (edge_seq, vertex_sig) in enumerate(cycles_list, 1):
            out.write(
                describe_graph(
                    edge_seq,
                    vertex_sig,
                    graph_num=i,
                    total=len(cycles_list),
                    compact=args.compact,
                )
            )
            out.write("\n")
    finally:
        if args.output:
            out.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
