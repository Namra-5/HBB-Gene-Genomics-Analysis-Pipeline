#!/usr/bin/env python3
"""
Script 01: FASTA Sequence Statistics
HBB Genomics Project — Task 1 Bonus Analysis

Parses the multi-sequence FASTA file and reports per-sequence statistics:
length, GC content, nucleotide composition, and a summary table.

Usage:
    python 01_fasta_stats.py
    python 01_fasta_stats.py --fasta path/to/Task1_hbb_seq.fasta
"""

import sys
import os
from collections import Counter

# ── Allow running from any directory ──────────────────────────────────────────
DEFAULT_FASTA = os.path.join(os.path.dirname(__file__), "..", "Task1_GenomicData", "Task1_hbb_seq.fasta")

# ── Population metadata ────────────────────────────────────────────────────────
POPULATION_MAP = {
    "NC_000011.10": ("GRCh38 Reference", "chromosome 11"),
    "MK476281.1":   ("Nzime",            "Cameroon"),
    "MK475906.1":   ("Bakoya",           "Cameroon"),
    "MK476491.1":   ("Yoruba",           "Nigeria"),
    "MK476080.1":   ("Bakota",           "Gabon"),
    "MK476374.1":   ("Kiga",             "Uganda"),
    "MK476446.1":   ("Mandenka",         "Senegal"),
}


def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    sequences = []
    current_header = None
    current_seq = []

    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq).upper()))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, "".join(current_seq).upper()))

    return sequences


def seq_stats(sequence):
    """Return a dict of statistics for a nucleotide sequence."""
    counts = Counter(sequence)
    length = len(sequence)
    gc = counts["G"] + counts["C"]
    at = counts["A"] + counts["T"]
    n_count = length - gc - at  # ambiguous bases

    return {
        "length": length,
        "A": counts["A"],
        "T": counts["T"],
        "G": counts["G"],
        "C": counts["C"],
        "N": n_count,
        "GC_count": gc,
        "GC_pct": round(gc / length * 100, 2),
        "AT_pct": round(at / length * 100, 2),
    }


def accession_from_header(header):
    """Extract the accession number from a FASTA header."""
    return header.split()[0].split(":")[0]


def main():
    # ── Locate FASTA file ──────────────────────────────────────────────────────
    fasta_path = sys.argv[2] if len(sys.argv) > 2 and sys.argv[1] == "--fasta" else DEFAULT_FASTA
    if not os.path.exists(fasta_path):
        # Try current directory fallback
        fasta_path = "Task1_hbb_seq.fasta"
    if not os.path.exists(fasta_path):
        print(f"[ERROR] Cannot find FASTA file. Please specify: python 01_fasta_stats.py --fasta <path>")
        sys.exit(1)

    print(f"\n{'='*70}")
    print("  HBB Gene FASTA Sequence Statistics — Task 1 Bonus Analysis")
    print(f"{'='*70}")
    print(f"  File: {fasta_path}\n")

    sequences = parse_fasta(fasta_path)
    print(f"  Sequences found: {len(sequences)}\n")

    all_stats = []
    for header, seq in sequences:
        acc = accession_from_header(header)
        pop, region = POPULATION_MAP.get(acc, ("Unknown", "—"))
        stats = seq_stats(seq)
        stats["accession"] = acc
        stats["population"] = pop
        stats["region"] = region
        all_stats.append(stats)

    # ── Per-sequence table ─────────────────────────────────────────────────────
    print(f"  {'Accession':<25} {'Population':<20} {'Region':<12} {'Length':>7} {'GC%':>6} {'A':>6} {'T':>6} {'G':>6} {'C':>6}")
    print(f"  {'-'*25} {'-'*20} {'-'*12} {'-'*7} {'-'*6} {'-'*6} {'-'*6} {'-'*6} {'-'*6}")

    total_len = 0
    for s in all_stats:
        print(f"  {s['accession']:<25} {s['population']:<20} {s['region']:<12} "
              f"{s['length']:>7,} {s['GC_pct']:>6.2f} "
              f"{s['A']:>6,} {s['T']:>6,} {s['G']:>6,} {s['C']:>6,}")
        total_len += s["length"]

    print(f"\n  Total nucleotides across all sequences: {total_len:,} bp")

    # ── GC content summary ─────────────────────────────────────────────────────
    gc_values = [s["GC_pct"] for s in all_stats]
    print(f"\n  GC Content Summary:")
    print(f"    Mean GC%  : {sum(gc_values)/len(gc_values):.2f}%")
    print(f"    Min GC%   : {min(gc_values):.2f}%  ({all_stats[gc_values.index(min(gc_values))]['accession']})")
    print(f"    Max GC%   : {max(gc_values):.2f}%  ({all_stats[gc_values.index(max(gc_values))]['accession']})")

    # ── Conservation check ─────────────────────────────────────────────────────
    pop_seqs = [(s["accession"], s) for s in all_stats if s["accession"] != "NC_000011.10"]
    ref = next(s for s in all_stats if s["accession"] == "NC_000011.10")
    all_identical_gc = all(s["GC_pct"] == ref["GC_pct"] for _, s in pop_seqs)

    print(f"\n  Conservation Check:")
    print(f"    Reference (NC_000011.10) GC%: {ref['GC_pct']}%")
    if all_identical_gc:
        print("    ✓ All population sequences share identical GC% with the reference")
    else:
        for acc, s in pop_seqs:
            diff = s["GC_pct"] - ref["GC_pct"]
            print(f"    {acc}: ΔGC% = {diff:+.2f}%")

    print(f"\n{'='*70}\n")

    # ── Try to generate a matplotlib figure (optional) ─────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        accessions = [s["accession"].split(":")[0] for s in all_stats]
        labels = [f"{s['population']}\n({s['accession'].split(':')[0]})" for s in all_stats]
        lengths = [s["length"] for s in all_stats]
        gc_pcts = [s["GC_pct"] for s in all_stats]
        colors = ["#e63946" if s["accession"].startswith("NC_") else "#457b9d" for s in all_stats]

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle("HBB Gene Sequence Statistics — 7 Sequences", fontsize=14, fontweight="bold")

        # Sequence lengths
        bars1 = axes[0].bar(range(len(labels)), lengths, color=colors, edgecolor="white", linewidth=0.5)
        axes[0].set_xticks(range(len(labels)))
        axes[0].set_xticklabels(labels, fontsize=7, rotation=20, ha="right")
        axes[0].set_ylabel("Sequence Length (bp)")
        axes[0].set_title("Sequence Lengths")
        axes[0].axhline(1608, color="#e63946", linestyle="--", linewidth=1, label="Reference (1,608 bp)")
        axes[0].axhline(1824, color="#457b9d", linestyle="--", linewidth=1, label="Population (1,824 bp)")
        axes[0].legend(fontsize=8)
        axes[0].set_ylim(1400, 2000)
        for bar, val in zip(bars1, lengths):
            axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10, f"{val:,}", ha="center", fontsize=7)

        # GC content
        bars2 = axes[1].bar(range(len(labels)), gc_pcts, color=colors, edgecolor="white", linewidth=0.5)
        axes[1].set_xticks(range(len(labels)))
        axes[1].set_xticklabels(labels, fontsize=7, rotation=20, ha="right")
        axes[1].set_ylabel("GC Content (%)")
        axes[1].set_title("GC Content per Sequence")
        axes[1].set_ylim(35, 50)
        for bar, val in zip(bars2, gc_pcts):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, f"{val:.1f}%", ha="center", fontsize=7)

        ref_patch = mpatches.Patch(color="#e63946", label="Reference Genome")
        pop_patch = mpatches.Patch(color="#457b9d", label="Population Vouchers")
        fig.legend(handles=[ref_patch, pop_patch], loc="lower center", ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.02))

        plt.tight_layout(rect=[0, 0.04, 1, 1])
        out_path = os.path.join(os.path.dirname(__file__), "..","Task7_Comparative Sequence Analysis and Pipeline Reproducibility", "figures", "sequence_statistics.png")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"  Figure saved to: {out_path}")
        plt.close()
    except ImportError:
        print("  [INFO] matplotlib not installed — skipping figure generation")
    except Exception as e:
        print(f"  [WARNING] Figure generation failed: {e}")


if __name__ == "__main__":
    main()
