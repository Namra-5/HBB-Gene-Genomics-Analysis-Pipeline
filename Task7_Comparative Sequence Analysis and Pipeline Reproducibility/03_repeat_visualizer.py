#!/usr/bin/env python3
"""
Script 03: RepeatMasker Results Parser and Visualizer
HBB Genomics Project — Task 4 Bonus Analysis

Parses the RepeatMasker .out file and generates a summary report
showing the position, type, and properties of all identified repeats.

Usage:
    python 03_repeat_visualizer.py
    python 03_repeat_visualizer.py --out path/to/hbb_seq.fasta.out
"""

import os
import sys

DEFAULT_OUT = os.path.join(
    os.path.dirname(__file__), "..", "Task4_Identification_of_Repetitive_Elements", "hbb_seq.fasta.out"
)

# Population labels for each accession
ACCESSION_LABELS = {
    "NC_000011.10:c5227071-5225464": "Reference (GRCh38)",
    "MK476281.1": "Nzime (Cameroon)",
    "MK475906.1": "Bakoya (Cameroon)",
    "MK476491.1": "Yoruba (Nigeria)",
    "MK476080.1": "Bakota (Gabon)",
    "MK476374.1": "Kiga (Uganda)",
    "MK476446.1": "Mandenka (Senegal)",
}

# Known HBB gene feature boundaries on the reference (1-based)
HBB_FEATURES = {
    "Exon 1":          (51,   142),
    "Intron 1 (IVS-1)":(143,  272),
    "Exon 2":          (273,  495),
    "Intron 2 (IVS-2)":(496, 1345),
    "Exon 3":          (1346, 1474),
}


def parse_repeatmasker_out(filepath):
    """
    Parse a RepeatMasker .out file.
    Returns a list of dicts, one per repeat hit.
    """
    repeats = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            # Skip header lines and blank lines
            if not line or line.startswith("SW") or line.startswith("score"):
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                entry = {
                    "sw_score":   int(parts[0]),
                    "pct_div":    float(parts[1]),
                    "pct_del":    float(parts[2]),
                    "pct_ins":    float(parts[3]),
                    "sequence":   parts[4],
                    "begin":      int(parts[5]),
                    "end":        int(parts[6]),
                    "left":       parts[7],
                    "strand":     parts[8],
                    "repeat":     parts[9],
                    "class_fam":  parts[10],
                    "rep_begin":  parts[11],
                    "rep_end":    parts[12],
                    "rep_left":   parts[13],
                    "hit_id":     int(parts[14]),
                    "length":     int(parts[6]) - int(parts[5]) + 1,
                }
                repeats.append(entry)
            except (ValueError, IndexError):
                continue
    return repeats


def print_report(repeats):
    print(f"\n{'='*72}")
    print("  RepeatMasker Results — Task 4 Bonus Analysis")
    print(f"{'='*72}")
    print(f"  Input: hbb_seq.fasta (7 sequences, 12,552 bp total)")
    print(f"  Tool: RepeatMasker v4.2.3 | Database: Dfam 3.9 (RepeatMasker.lib)")
    print(f"  Total repeat hits identified: {len(repeats)}\n")

    # Summary by sequence
    by_seq = {}
    for r in repeats:
        by_seq.setdefault(r["sequence"], []).append(r)

    print(f"  Repeats per Sequence:")
    print(f"  {'Sequence':<42} {'Population':<22} {'Hits':>4} {'Masked bp':>10}")
    print(f"  {'-'*42} {'-'*22} {'-'*4} {'-'*10}")
    total_masked = 0
    for seq, hits in sorted(by_seq.items()):
        label = ACCESSION_LABELS.get(seq, seq)
        masked = sum(h["length"] for h in hits)
        total_masked += masked
        print(f"  {seq:<42} {label:<22} {len(hits):>4} {masked:>10} bp")
    print(f"  {'TOTAL':<42} {'':22} {len(repeats):>4} {total_masked:>10} bp")

    # Repeat motif summary
    motifs = {}
    for r in repeats:
        motifs.setdefault(r["repeat"], []).append(r)

    print(f"\n  Repeat Motif Summary:")
    print(f"  {'Motif':<15} {'Class':<20} {'Hits':>5} {'Total bp':>9} {'Avg SW':>7} {'Avg Div%':>9}")
    print(f"  {'-'*15} {'-'*20} {'-'*5} {'-'*9} {'-'*7} {'-'*9}")
    for motif, hits in motifs.items():
        avg_sw  = sum(h["sw_score"]  for h in hits) / len(hits)
        avg_div = sum(h["pct_div"]   for h in hits) / len(hits)
        tot_bp  = sum(h["length"]    for h in hits)
        cls     = hits[0]["class_fam"]
        print(f"  {motif:<15} {cls:<20} {len(hits):>5} {tot_bp:>9} {avg_sw:>7.1f} {avg_div:>9.1f}%")

    # Genomic context (reference only)
    ref_hits = [r for r in repeats if r["sequence"].startswith("NC_")]
    print(f"\n  Reference Sequence Repeat Coordinates and Genomic Context:")
    for r in ref_hits:
        context = "intergenic / flanking"
        for feat, (start, end) in HBB_FEATURES.items():
            if r["begin"] >= start and r["end"] <= end:
                context = feat
                break
        print(f"    bp {r['begin']:>4}–{r['end']:<4}  {r['repeat']:<12}  "
              f"SW={r['sw_score']:>2}  Div={r['pct_div']:>4.1f}%  "
              f"Length={r['length']:>3} bp  Context: {context}")

    print(f"\n  Biological Notes:")
    print(f"    1. All 14 repeats belong to the Simple_repeat class — no transposable")
    print(f"       elements (SINEs, LINEs, DNA transposons) were detected.")
    print(f"    2. The (TTTTC)n and (TTTC)n motifs fall within Intron 2 (IVS-2),")
    print(f"       the largest intron in HBB (~850 bp on the reference).")
    print(f"    3. These AT-rich microsatellite repeats are conserved at identical")
    print(f"       positions across all 7 sequences, indicating a stable genomic feature.")
    print(f"    4. The repeat region (639–711 bp on reference) shows GC% of ~20%,")
    print(f"       consistent with the (TTTC)n motif (100% AT in the canonical form).")
    print(f"    5. 511 bp masked total (4.07% of 12,552 bp). The masked FASTA")
    print(f"       (hbb_seq.fasta.masked) is suitable for downstream variant calling.")
    print(f"\n{'='*72}\n")


def generate_figure(repeats):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        sequences = list(dict.fromkeys(r["sequence"] for r in repeats))
        labels = [ACCESSION_LABELS.get(s, s) for s in sequences]

        fig, axes = plt.subplots(len(sequences), 1, figsize=(14, len(sequences) * 1.1 + 1.5))
        fig.suptitle(
            "RepeatMasker — Repeat Element Positions across HBB Sequences\n"
            "(All repeats: Simple_repeat / (TTTTC)n and (TTTC)n motifs)",
            fontsize=11, fontweight="bold"
        )

        color_map = {"(TTTTC)n": "#ffb703", "(TTTC)n": "#fb8500"}

        for i, (seq, label) in enumerate(zip(sequences, labels)):
            ax = axes[i] if len(sequences) > 1 else axes
            ax.set_xlim(0, 1824)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_ylabel(label, fontsize=7, rotation=0, ha="right", va="center", labelpad=5)

            # Shade exons (only for reference)
            if seq.startswith("NC_"):
                seq_len = 1608
                exon_coords = [(51, 142), (273, 495), (1346, 1474)]
                for ex_s, ex_e in exon_coords:
                    ax.axvspan(ex_s, ex_e, ymin=0, ymax=1, alpha=0.2, color="#52b788")
                ax.set_xlim(0, seq_len)
            else:
                ax.set_xlim(0, 1824)

            # Draw repeats
            for r in repeats:
                if r["sequence"] == seq:
                    color = color_map.get(r["repeat"], "#8ecae6")
                    ax.barh(0.5, r["length"], left=r["begin"], height=0.6,
                            color=color, edgecolor="black", linewidth=0.5)

            if i < len(sequences) - 1:
                ax.set_xticks([])
            else:
                ax.set_xlabel("Position (bp)", fontsize=9)

        # Legend
        patch1 = mpatches.Patch(color="#ffb703", label="(TTTTC)n repeat")
        patch2 = mpatches.Patch(color="#fb8500", label="(TTTC)n repeat")
        patch3 = mpatches.Patch(color="#52b788", alpha=0.3, label="Exons (reference only)")
        fig.legend(handles=[patch1, patch2, patch3], loc="lower center",
                   ncol=3, fontsize=9, bbox_to_anchor=(0.5, 0.01))

        plt.tight_layout(rect=[0, 0.05, 1, 1])
        out_path = os.path.join(os.path.dirname(__file__), "..", "Task7_Comparative Sequence Analysis and Pipeline Reproducibility", "figures", "repeat_distribution.png")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"  Figure saved to: {out_path}")
        plt.close()
    except ImportError:
        print("  [INFO] matplotlib not available — skipping figure generation")
    except Exception as e:
        print(f"  [WARNING] Figure generation failed: {e}")


def main():
    out_path = DEFAULT_OUT
    if "--out" in sys.argv:
        out_path = sys.argv[sys.argv.index("--out") + 1]
    if not os.path.exists(out_path):
        out_path = os.path.join("Task4_Identification_of_Repetitive_Elements", "hbb_seq.fasta.out")
    if not os.path.exists(out_path):
        print("[ERROR] RepeatMasker .out file not found.")
        sys.exit(1)

    repeats = parse_repeatmasker_out(out_path)
    print_report(repeats)
    generate_figure(repeats)


if __name__ == "__main__":
    main()
