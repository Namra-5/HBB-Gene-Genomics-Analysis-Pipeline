#!/usr/bin/env python3
"""
Script 02: GC Content Sliding Window Analysis
HBB Genomics Project — Task 3 & 4 Bonus Analysis

Computes GC content in a sliding window across the HBB reference sequence
and overlays the positions of RepeatMasker-identified simple repeats.

This analysis reveals:
  - GC-rich regions (exons, regulatory elements)
  - AT-rich regions (introns, repeat regions)
  - The exact location of (TTTC)n microsatellite clusters

Usage:
    python 02_gc_content_analysis.py
    python 02_gc_content_analysis.py --window 50 --step 10
"""

import os
import sys

# ── Defaults ───────────────────────────────────────────────────────────────────
DEFAULT_FASTA = os.path.join(os.path.dirname(__file__), "..", "Task1_GenomicData", "Task1_hbb_seq.fasta")
WINDOW_SIZE   = 50   # bp sliding window
STEP_SIZE     = 10   # bp step

# ── Known HBB gene features (positions on the 1,608 bp reference) ─────────────
HBB_FEATURES = {
    "Exon 1":  (51,  142),    # encodes aa 1–30, includes sickle-cell codon 6
    "Intron 1 (IVS-1)": (143, 272),
    "Exon 2":  (273, 495),    # encodes aa 31–104, haem-binding region
    "Intron 2 (IVS-2)": (496, 1345),
    "Exon 3":  (1346, 1474),  # encodes aa 105–146
}

# ── RepeatMasker results (from .out file) ──────────────────────────────────────
REPEATS = [
    # (start, end, motif, SW_score, divergence)
    (639, 674, "(TTTTC)n", 22, 29.6),
    (675, 711, "(TTTC)n",  14, 12.5),
]


def parse_fasta(filepath):
    sequences = {}
    current_header = None
    current_seq = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq).upper()
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = "".join(current_seq).upper()
    return sequences


def sliding_gc(sequence, window=50, step=10):
    """Compute GC% in a sliding window. Returns (midpoints, gc_values)."""
    midpoints, gc_values = [], []
    for start in range(0, len(sequence) - window + 1, step):
        window_seq = sequence[start:start + window]
        gc = (window_seq.count("G") + window_seq.count("C")) / window * 100
        midpoints.append(start + window // 2)
        gc_values.append(gc)
    return midpoints, gc_values


def print_text_report(sequence, midpoints, gc_values):
    """Print a text-based GC content summary."""
    print(f"\n{'='*70}")
    print("  HBB Reference Sequence — GC Content Analysis")
    print(f"{'='*70}")
    print(f"  Sequence: NC_000011.10 (GRCh38 HBB locus, 1,608 bp)")
    print(f"  Window size: {WINDOW_SIZE} bp | Step: {STEP_SIZE} bp")
    print(f"  Windows computed: {len(midpoints)}\n")

    overall_gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    print(f"  Overall GC%   : {overall_gc:.2f}%")
    print(f"  Mean window GC: {sum(gc_values)/len(gc_values):.2f}%")
    print(f"  Min GC%       : {min(gc_values):.1f}%  (around position {midpoints[gc_values.index(min(gc_values))]})")
    print(f"  Max GC%       : {max(gc_values):.1f}%  (around position {midpoints[gc_values.index(max(gc_values))]})")

    print(f"\n  Known HBB Gene Features:")
    for feature, (start, end) in HBB_FEATURES.items():
        feature_seq = sequence[start-1:end]
        gc = (feature_seq.count("G") + feature_seq.count("C")) / len(feature_seq) * 100
        print(f"    {feature:<25} bp {start:>4}–{end:<4}  ({end-start+1:>4} bp)  GC%: {gc:.1f}%")

    print(f"\n  RepeatMasker Simple Repeats (reference coordinates):")
    for start, end, motif, sw, div in REPEATS:
        rpt_seq = sequence[start-1:end]
        gc = (rpt_seq.count("G") + rpt_seq.count("C")) / len(rpt_seq) * 100
        print(f"    {motif:<12}  bp {start:>3}–{end:<3}  SW={sw}  Div={div}%  GC%: {gc:.1f}%  ← AT-rich repeat region")

    print(f"\n  Biological Interpretation:")
    print("    • Exons 2 and 3 (haem-binding and C-terminal regions) show higher GC%,")
    print("      consistent with codon composition constraints in functional protein-coding regions.")
    print("    • The (TTTC)n repeat region at ~639–711 bp falls within Intron 2 (IVS-2),")
    print("      in an AT-rich environment typical of intronic simple sequence repeats.")
    print("    • AT-rich introns are consistent with nucleosome-depleted, accessible chromatin")
    print("      in the HBB locus, which is regulated by the upstream Locus Control Region (LCR).")
    print(f"\n{'='*70}\n")


def main():
    # ── Parse arguments ────────────────────────────────────────────────────────
    global WINDOW_SIZE, STEP_SIZE
    if "--window" in sys.argv:
        WINDOW_SIZE = int(sys.argv[sys.argv.index("--window") + 1])
    if "--step" in sys.argv:
        STEP_SIZE = int(sys.argv[sys.argv.index("--step") + 1])

    fasta_path = DEFAULT_FASTA
    if not os.path.exists(fasta_path):
        fasta_path = "Task1_hbb_seq.fasta"
    if not os.path.exists(fasta_path):
        print("[ERROR] FASTA file not found.")
        sys.exit(1)

    sequences = parse_fasta(fasta_path)

    # Use the reference sequence (NC_000011.10)
    ref_key = next((k for k in sequences if k.startswith("NC_000011.10")), None)
    if not ref_key:
        ref_key = list(sequences.keys())[0]
    ref_seq = sequences[ref_key]

    midpoints, gc_values = sliding_gc(ref_seq, WINDOW_SIZE, STEP_SIZE)
    print_text_report(ref_seq, midpoints, gc_values)

    # ── Matplotlib figure ──────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        fig, ax = plt.subplots(figsize=(14, 5))

        # GC line
        ax.plot(midpoints, gc_values, color="#1d3557", linewidth=1.2, alpha=0.85, label=f"GC% (window={WINDOW_SIZE} bp)")
        overall_gc = (ref_seq.count("G") + ref_seq.count("C")) / len(ref_seq) * 100
        ax.axhline(overall_gc, color="#e63946", linestyle="--", linewidth=1, alpha=0.7, label=f"Overall GC% ({overall_gc:.1f}%)")

        # Shade exons
        exon_colors = ["#a8dadc", "#52b788", "#95d5b2"]
        for i, (feature, (start, end)) in enumerate(HBB_FEATURES.items()):
            if "Exon" in feature:
                color = exon_colors[i // 2]
                ax.axvspan(start, end, alpha=0.25, color=color, label=feature)

        # Mark repeats
        for start, end, motif, sw, div in REPEATS:
            ax.axvspan(start, end, alpha=0.4, color="#ffb703", label=f"Repeat: {motif}")

        ax.set_xlabel("Position in HBB reference sequence (bp)", fontsize=11)
        ax.set_ylabel("GC Content (%)", fontsize=11)
        ax.set_title(f"HBB Locus GC Content — Sliding Window Analysis\n(NC_000011.10, GRCh38, window={WINDOW_SIZE} bp, step={STEP_SIZE} bp)", fontsize=12)
        ax.set_xlim(0, len(ref_seq))
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3, linestyle="--")

        # De-duplicate legend labels
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=8, ncol=2)

        plt.tight_layout()
        out_path = os.path.join(os.path.dirname(__file__), "..", "Task7_Comparative Sequence Analysis and Pipeline Reproducibility","figures", "gc_content_sliding_window.png")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"  Figure saved to: {out_path}")
        plt.close()

    except ImportError:
        print("  [INFO] matplotlib not available — skipping figure generation")
    except Exception as e:
        print(f"  [WARNING] Figure generation failed: {e}")


if __name__ == "__main__":
    main()
