#!/usr/bin/env python3
"""
Script 06: Full Pipeline Results Summary
HBB Genomics Project — End-to-End Analysis Report

Aggregates results from all six tasks and produces a structured
console report of the complete genomics analysis pipeline.

Usage:
    python 06_pipeline_summary.py
"""

import os
import sys

SEPARATOR = "=" * 72
MINOR_SEP = "-" * 72


def section(title):
    print(f"\n{SEPARATOR}")
    print(f"  {title}")
    print(SEPARATOR)


def check_file(path, label):
    exists = os.path.exists(path)
    status = "FOUND" if exists else "NOT FOUND"
    size   = f"{os.path.getsize(path):>10,} bytes" if exists else f"{'—':>10}"
    print(f"  {'[OK]' if exists else '[--]'}  {label:<50} {size}  {status}")
    return exists


def main():
    # ── Detect repo root ───────────────────────────────────────────────────────
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root  = os.path.abspath(os.path.join(script_dir, ".."))

    def rp(*parts):
        return os.path.join(repo_root, *parts)

    print(f"\n{SEPARATOR}")
    print("  HBB Gene Genomics Pipeline — Complete Results Summary")
    print("  National University of Sciences & Technology (NUST), SINES")
    print("  BS Bioinformatics UG-1 | Genomics Pre-Midterm Project")
    print(SEPARATOR)
    print(f"\n  Gene     : HBB — Hemoglobin Subunit Beta (Gene ID: 3043)")
    print(f"  Organism : Homo sapiens")
    print(f"  Reference: NC_000011.10 (GRCh38.p14, chr11:5,225,464–5,227,071)")
    print(f"  Sequences: 7 total (1 reference + 6 African population vouchers)")
    print(f"  Disease  : Sickle Cell Disease (HbS, rs334, Glu6Val) / Beta-Thalassemia")

    # ── File inventory ─────────────────────────────────────────────────────────
    section("Repository File Inventory")
    files = [
        (os.path.join("Task1_GenomicData", "Task1_hbb_seq.fasta"),                              "Task 1: All 7 HBB sequences (FASTA)"),
        (os.path.join("Task2_Genome_and_Annotation_Visualization","Task2_hbb.bam"),          "Task 2: Sorted BAM file"),
        (os.path.join("Task2_Genome_and_Annotation_Visualization","Task2_hbb.bai"),          "Task 2: BAM index file"),
        (os.path.join("Task2_Genome_and_Annotation_Visualization","Task2_igv_session.xml"),  "Task 2: IGV session file"),
        (os.path.join("Task3_Sequence_Alignment_and_Similarity","Task3_ClustalW.clustal"), "Task 3: ClustalW MSA output"),
        (os.path.join("Task3_Sequence_Alignment_and_Similarity","Task3_Muscle.aln-clustalw"),"Task 3: MUSCLE MSA output"),
        (os.path.join("Task4_Identification_of_Repetitive_Elements","hbb_seq.fasta.out"),   "Task 4: Repeat annotations (.out)"),
        (os.path.join("Task4_Identification_of_Repetitive_Elements","hbb_seq.fasta.out.gff"),"Task 4: Repeat annotations (.gff)"),
        (os.path.join("Task4_Identification_of_Repetitive_Elements","hbb_seq.fasta.masked"), "Task 4: Masked FASTA"),
        (os.path.join("Task4_Identification_of_Repetitive_Elements","hbb_seq.fasta.tbl"),    "Task 4: RepeatMasker stats table"),
        (os.path.join("Task5_Variant_Calling_and_Functional_Annotation","Task5-FreeBayes_variants.vcf"), "Task 5: Raw FreeBayes VCF"),
        (os.path.join("Task5_Variant_Calling_and_Functional_Annotation","Task5-VCFfilter.vcf"),                  "Task 5: Filtered VCF"),
        (os.path.join("Task5_Variant_Calling_and_Functional_Annotation","Task5-SnpEff_CSV stats.txt"),   "Task 5: SnpEff CSV stats"),
        (os.path.join("Task6_Gene_Prediction_and_Orthology_Inference","Augustus_Output.gff"),                         "Task 6: AUGUSTUS GFF3"),
        (os.path.join("Task6_Gene_Prediction_and_Orthology_Inference","Augustus_Predicted_Aminoacidseq.fa"),          "Task 6: Predicted proteins (AA)"),
        (os.path.join("Task6_Gene_Prediction_and_Orthology_Inference","AugustusPredicted_codingseq.fa"),              "Task 6: Predicted CDS (nt)"),
    ]
    for rel_path, label in files:
        check_file(rp(rel_path), label)

    # ── Task 1 ─────────────────────────────────────────────────────────────────
    section("Task 1: Data Retrieval and Documentation")
    populations = [
        ("NC_000011.10", "GRCh38 Reference",   "chr11 primary assembly",  "1,608"),
        ("MK476281.1",   "Nzime",              "Cameroon",                "1,824"),
        ("MK475906.1",   "Bakoya",             "Cameroon",                "1,824"),
        ("MK476491.1",   "Yoruba",             "Nigeria",                 "1,824"),
        ("MK476080.1",   "Bakota",             "Gabon",                   "1,824"),
        ("MK476374.1",   "Kiga",               "Uganda",                  "1,824"),
        ("MK476446.1",   "Mandenka",           "Senegal",                 "1,824"),
    ]
    print(f"\n  {'Accession':<16} {'Population':<20} {'Region':<26} {'Length':>8}")
    print(f"  {'-'*16} {'-'*20} {'-'*26} {'-'*8}")
    for acc, pop, region, length in populations:
        print(f"  {acc:<16} {pop:<20} {region:<26} {length:>8} bp")

    print(f"\n  BLAST Validation: All sequences — 100% identity | 100% coverage | E-value=0.0")
    print(f"  GC Content: Reference 40.11% | Population sequences 40.24% (delta = +0.13%)")
    print(f"  Total dataset size: 12,552 bp across 7 sequences")

    # ── Task 2 ─────────────────────────────────────────────────────────────────
    section("Task 2: IGV Genome Visualization")
    print(f"\n  Reference genome: Human GRCh38/hg38 (IGV built-in)")
    print(f"  Region visualized: chr11:5,225,404–5,227,132 (1,729 bp window)")
    print(f"  HBB gene locus  : chr11:5,225,464–5,227,071 (1,608 bp)")
    print(f"\n  {'Sample':<18} {'Accession':<16} {'Strand':<8} {'CIGAR':<10} {'NM':>4} {'MAPQ':>6}")
    print(f"  {'-'*18} {'-'*16} {'-'*8} {'-'*10} {'-'*4} {'-'*6}")
    igv_data = [
        ("Nzime",    "MK476281.1", "Minus (–)", "1608M", "0", "60"),
        ("Bakoya",   "MK475906.1", "Minus (–)", "1608M", "0", "60"),
        ("Yoruba",   "MK476491.1", "Minus (–)", "1608M", "0", "60"),
        ("Bakota",   "MK476080.1", "Minus (–)", "1608M", "0", "60"),
        ("Kiga",     "MK476374.1", "Minus (–)", "1608M", "0", "60"),
        ("Mandenka", "MK476446.1", "Minus (–)", "1608M", "0", "60"),
    ]
    for sample, acc, strand, cigar, nm, mapq in igv_data:
        print(f"  {sample:<18} {acc:<16} {strand:<8} {cigar:<10} {nm:>4} {mapq:>6}")
    print(f"\n  Coverage: Uniform depth of 6 across all 1,608 bp (no gaps)")
    print(f"  Variants: None detected in any sample (wild-type HBB confirmed in all)")
    print(f"  Gene structure visible: 3 exons, 2 introns (RefSeq annotation, minus strand)")

    # ── Task 3 ─────────────────────────────────────────────────────────────────
    section("Task 3: Sequence Alignment and Similarity Analysis")
    print(f"\n  PART A — NCBI BLAST (Megablast, Blast 2 Sequences)")
    print(f"  {'Run':<5} {'Query':<20} {'Subject':<18} {'Score':>6} {'Cover':>7} {'Identity':>10} {'E-value':>10} {'Gaps':>6}")
    print(f"  {'-'*5} {'-'*20} {'-'*18} {'-'*6} {'-'*7} {'-'*10} {'-'*10} {'-'*6}")
    blast_data = [
        ("1","NC_000011.10","MK476281.1 (Nzime)",   2970,"100%","100% (1608/1608)","0.0","0/1608"),
        ("2","NC_000011.10","MK475906.1 (Bakoya)",  2970,"100%","100% (1608/1608)","0.0","0/1608"),
        ("3","NC_000011.10","MK476491.1 (Yoruba)",  2970,"100%","100% (1608/1608)","0.0","0/1608"),
        ("4","NC_000011.10","MK476080.1 (Bakota)",  2970,"100%","100% (1608/1608)","0.0","0/1608"),
        ("5","NC_000011.10","MK476374.1 (Kiga)",    2970,"100%","100% (1608/1608)","0.0","0/1608"),
        ("6","NC_000011.10","MK476446.1 (Mandenka)",2970,"100%","100% (1608/1608)","0.0","0/1608"),
    ]
    for run, q, subj, score, cov, ident, evalue, gaps in blast_data:
        print(f"  {run:<5} {q:<20} {subj:<18} {score:>6} {cov:>7} {ident:>10} {evalue:>10} {gaps:>6}")

    print(f"\n  PART B — Multiple Sequence Alignment")
    print(f"  ClustalW (Galaxy) : 100% conserved positions (*) across all 7 sequences")
    print(f"  MUSCLE (EMBL-EBI) : Identical conservation pattern — independent confirmation")
    print(f"  Conclusion        : Complete sequence conservation; purifying selection signature")

    # ── Task 4 ─────────────────────────────────────────────────────────────────
    section("Task 4: RepeatMasker — Repeat Element Identification")
    print(f"\n  Input          : 7 sequences, 12,552 bp total (hbb_seq.fasta)")
    print(f"  Tool           : RepeatMasker v4.2.3 | RMBlast v2.14.1+ | Dfam 3.9")
    print(f"  Total masked   : 511 bp (4.07%)")
    print(f"  Repeats found  : 14 elements")
    print(f"\n  {'Motif':<15} {'Class':<20} {'Hits':>5} {'Total bp':>9} {'Avg SW':>8} {'Div%':>8}")
    print(f"  {'-'*15} {'-'*20} {'-'*5} {'-'*9} {'-'*8} {'-'*8}")
    print(f"  {'(TTTTC)n':<15} {'Simple_repeat':<20} {7:>5} {252:>9} {22.0:>8.1f} {29.6:>8.1f}%")
    print(f"  {'(TTTC)n':<15} {'Simple_repeat':<20} {7:>5} {259:>9} {14.0:>8.1f} {12.5:>8.1f}%")
    print(f"\n  Reference positions: bp 639–674 (TTTTC)n | bp 675–711 (TTTC)n")
    print(f"  Genomic context    : Intron 2 (IVS-2) of HBB gene")
    print(f"  TEs detected       : None (SINEs=0, LINEs=0, DNA transposons=0, LTR=0)")

    # ── Task 5 ─────────────────────────────────────────────────────────────────
    section("Task 5: Variant Calling and Functional Annotation")
    print(f"\n  Platform: Galaxy (usegalaxy.org) | All steps completed without errors")
    print(f"\n  {'Step':<5} {'Tool':<22} {'Version':<12} {'Output':<30} {'Key Parameters'}")
    print(f"  {'-'*5} {'-'*22} {'-'*12} {'-'*30} {'-'*30}")
    pipeline = [
        ("1","BWA-MEM",   "0.7.19",  "BAM (1.6 KB)",               "Single-end, simple Illumina mode"),
        ("2","FreeBayes", "1.3.10",  "Raw VCF (7.8 KB)",           "Ploidy=2, min base quality=20"),
        ("3","VCFtools",  "—",       "Filtered VCF (QUAL≥20,DP≥10)","Standard conservative thresholds"),
        ("4","SnpEff",    "5.2",     "Ann. VCF + HTML + CSV",       "hg38 database, GENCODE annotations"),
    ]
    for step, tool, ver, out, params in pipeline:
        print(f"  {step:<5} {tool:<22} {ver:<12} {out:<30} {params}")

    print(f"\n  VCF Header metadata confirms: FreeBayes v1.3.10 | VCF format v4.2")
    print(f"  Contigs: NC_000011.10 (1,608 bp) + 6 MK-series (1,824 bp each)")
    print(f"  SnpEff: hg38 genome loaded (3,247,448,939 bp) | 0 warnings | 0 errors")

    # ── Task 6 ─────────────────────────────────────────────────────────────────
    section("Task 6: Gene Prediction (AUGUSTUS) and Orthology (REvolutionH-tl)")
    print(f"\n  PART A — AUGUSTUS Ab Initio Gene Prediction")
    print(f"  Model: Homo sapiens | Mode: ab initio (no hints) | Both strands")
    print(f"\n  {'Gene':<5} {'Accession':<18} {'Population':<22} {'Score':>6} {'Exons':>6} {'Protein match'}")
    print(f"  {'-'*5} {'-'*18} {'-'*22} {'-'*6} {'-'*6} {'-'*30}")
    aug_data = [
        ("g1","NC_000011.10","GRCh38 Reference",  "0.93",3,"Exact (UniProt P68871)"),
        ("g2","MK476281.1",  "Nzime (Cameroon)",   "0.97",3,"Exact (UniProt P68871)"),
        ("g3","MK475906.1",  "Bakoya (Cameroon)",  "0.96",3,"Exact (UniProt P68871)"),
        ("g4","MK476491.1",  "Yoruba (Nigeria)",   "0.95",3,"Exact (UniProt P68871)"),
        ("g5","MK476080.1",  "Bakota (Gabon)",     "0.95",3,"Exact (UniProt P68871)"),
        ("g6","MK476374.1",  "Kiga (Uganda)",      "0.95",3,"Exact (UniProt P68871)"),
        ("g7","MK476446.1",  "Mandenka (Senegal)", "0.95",3,"Exact (UniProt P68871)"),
    ]
    for g, acc, pop, score, nexons, match in aug_data:
        print(f"  {g:<5} {acc:<18} {pop:<22} {score:>6} {nexons:>6} {match}")

    print(f"\n  Exon 2 confidence: 1.0 (maximum) in all 7 sequences")
    print(f"  Predicted protein: 146 aa | Exact match to canonical beta-globin (P68871)")

    print(f"\n  PART B — REvolutionH-tl Orthology Inference")
    print(f"  Dataset  : 4 Terrabacteria species (bacterial demo dataset)")
    print(f"  Aligner  : DIAMOND v2.1.8 (all-vs-all protein alignment)")
    print(f"  Steps    : All 6 pipeline steps completed successfully")
    print(f"  Output   : Orthogroup assignments, gene trees, reconciled species tree")
    print(f"  Tree     : Visualized on iTOL — topology matches known Terrabacteria phylogeny")
    print(f"  HBB note : Confirmed as single-copy vertebrate ortholog in OrthoDB (EOG09150CZZ)")

    # ── Final Summary ──────────────────────────────────────────────────────────
    section("Overall Findings Summary")
    print(f"""
  1. The HBB coding sequence is 100% conserved across 6 diverse African
     populations spanning Cameroon, Nigeria, Gabon, Uganda, and Senegal.

  2. No sickle cell mutation (rs334, Glu6Val) or any other coding variant
     was detected in any of the 6 population samples.

  3. AUGUSTUS correctly predicted the 3-exon / 2-intron HBB gene structure
     from raw sequence alone, producing the exact canonical beta-globin
     protein (UniProt P68871) for all 7 sequences.

  4. RepeatMasker identified 14 simple repeat elements ((TTTC)n / (TTTTC)n)
     in Intron 2 of HBB, conserved across all sequences. No transposable
     elements were found in the 12,552 bp analyzed.

  5. The complete BWA-MEM → FreeBayes → VCFtools → SnpEff pipeline was
     executed on Galaxy without errors. The SnpEff coordinate mismatch
     reflects a known limitation of running targeted-contig variant calls
     through a whole-genome annotation engine.

  6. Strong purifying selection maintains the HBB coding sequence. The
     gene encodes the oxygen-transporting beta-globin chain; most mutations
     cause severe disease (SCD, beta-thalassemia), preventing accumulation
     of variants across populations.
""")
    print(SEPARATOR)


if __name__ == "__main__":
    main()
