#!/usr/bin/env python3
"""
Read the Kelley penalty table that Schrödinger's conformer_cluster.py writes
(<name>_cluster_ligand1.txt) and print the number of clusters with the lowest penalty at 2 or more
clusters. Optionally save a plot of the penalty curve. Called by 02_cluster_conformers.sh.

This automates the earlier manual fix ("get_second_lowest_Kelley_penalty_MMclustering.py"): when the
1-cluster row has the lowest penalty, use the lowest penalty among the remaining rows. Taking the minimum
over k >= 2 gives the same answer whether or not k = 1 is the overall minimum.

    python kelley_schrodinger.py mol233_cluster_ligand1.txt                     -> prints e.g. 5
    python kelley_schrodinger.py mol233_cluster_ligand1.txt --plot mol233.png --title mol233

The penalty column is found by its header "Kelley Penalty"; if that header is missing, the 6th column
is used, as in the old script.
"""
import argparse, csv, sys


def read_table(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r]
    head = [h.strip() for h in rows[0]]
    col = next((i for i, h in enumerate(head) if h.lower() == "kelley penalty"), 5)
    ks, pen = [], []
    for r in rows[1:]:
        try:
            k, p = int(float(r[0])), float(r[col])
        except (ValueError, IndexError):
            continue                                  # rows without a penalty value are skipped
        ks.append(k); pen.append(p)
    if not ks: raise ValueError(f"no Kelley penalty values found in {path}")
    return ks, pen


def best_k(ks, pen, min_k=2):
    cand = [(p, k) for k, p in zip(ks, pen) if k >= min_k]
    if not cand: raise ValueError("no rows with 2 or more clusters")
    return min(cand)[1]                               # ties -> fewer clusters


def plot(ks, pen, k_best, k_default, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt, matplotlib.ticker
    blue, orange, grey, muted, grid = "#2a78d6", "#eb6834", "#9b9a95", "#52514e", "#e6e5e0"
    fig, ax = plt.subplots(figsize=(5.2, 3), dpi=130)
    ax.grid(color=grid, lw=.7); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.plot(ks, pen, color=blue, marker="o", ms=3, lw=1.4)
    if k_default != k_best:                           # Schrödinger's own (default) choice, e.g. 1 cluster
        ax.plot([k_default], [pen[ks.index(k_default)]], marker="o", ms=7, mfc="white", mec=grey, mew=1.5, zorder=3)
    ax.plot([k_best], [pen[ks.index(k_best)]], marker="o", ms=7, color=orange, zorder=4)
    ax.axvline(k_best, color=orange, lw=1, ls="--")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_xlabel("Number of clusters", color=muted); ax.set_ylabel("Kelley penalty (Schrödinger)", color=muted)
    sub = f"used: {k_best}" + (f" (default minimum was {k_default}, open circle)" if k_default != k_best else " (default minimum)")
    ax.set_title(f"{title}\n{sub}", loc="left", fontsize=9, fontweight="bold")
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--plot", help="save the penalty curve to this .png")
    ap.add_argument("--title", default="")
    a = ap.parse_args()
    ks, pen = read_table(a.table)
    k = best_k(ks, pen)
    k_default = best_k(ks, pen, min_k=1)
    if a.plot:
        try:
            plot(ks, pen, k, k_default, a.plot, a.title)
        except ImportError:
            print("matplotlib not available; no plot written", file=sys.stderr)
    print(k)


if __name__ == "__main__":
    main()
