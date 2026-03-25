# HBB Gene Genomics Analysis Pipeline

**Genomics Pre-Midterm Project | BS Bioinformatics UG-1**  
School of Interdisciplinary Engineering and Sciences (SINES)  
National University of Sciences and Technology (NUST), Islamabad

| Field | Detail |
|-------|--------|
| Authors | Namra Basharat, Ghania Munir, Hania Fahad, Nawal Babar |
| Roll Numbers | 476203, 460673, 455132, 459508 |
| Gene | HBB — Hemoglobin Subunit Beta (Gene ID: 3043) |
| Organism | *Homo sapiens* |
| Reference | NC_000011.10 (GRCh38.p14, chr11:5,225,464–5,227,071) |
| Date | March 2026 |

---

## Project Overview

This repository contains the complete data, output files, and analysis scripts for a six-task genomics pipeline applied to the human **HBB (Hemoglobin Subunit Beta)** gene. The HBB gene was selected because a single point mutation at codon 6 (Glu6Val, rs334) causes Sickle Cell Disease (SCD), one of the most prevalent inherited disorders worldwide, with the highest burden in sub-Saharan Africa.

Seven sequences were analyzed: one GRCh38 reference sequence and six population-level voucher sequences from African populations in Cameroon, Nigeria, Gabon, Uganda, and Senegal. The analysis covers data retrieval, genome visualization, pairwise and multiple sequence alignment, repeat element identification, variant calling and annotation, gene structure prediction, and orthology inference.

---

## Dataset

| No. | Accession | Population | Country | Length |
|-----|-----------|-----------|---------|--------|
| 1 | NC_000011.10 | GRCh38 Reference | — | 1,608 bp |
| 2 | MK476281.1 | Nzime | Cameroon | 1,824 bp |
| 3 | MK475906.1 | Bakoya | Cameroon | 1,824 bp |
| 4 | MK476491.1 | Yoruba | Nigeria | 1,824 bp |
| 5 | MK476080.1 | Bakota | Gabon | 1,824 bp |
| 6 | MK476374.1 | Kiga | Uganda | 1,824 bp |
| 7 | MK476446.1 | Mandenka | Senegal | 1,824 bp |

All sequences retrieved from NCBI Nucleotide and validated by BLASTN (100% identity, 100% query coverage, E-value = 0.0).

---

## Analysis Pipeline

```
NCBI Nucleotide Database
        |
        v
Task 1: Data Retrieval
  NCBI search + BLASTN validation
  7 HBB sequences (FASTA)
        |
        v
Task 2: Genome Visualization
  IGV Desktop + hg38 reference
  BAM file prepared via Galaxy
        |
        v
Task 3: Sequence Alignment
  BLAST (pairwise) + ClustalW + MUSCLE (MSA)
        |
        v
Task 4: Repeat Identification
  RepeatMasker v4.2.3 + Dfam 3.9
  Masked FASTA for downstream analysis
        |
        v
Task 5: Variant Calling and Annotation
  BWA-MEM -> FreeBayes -> VCFtools -> SnpEff
  (Galaxy platform)
        |
        v
Task 6: Gene Prediction and Orthology
  AUGUSTUS (ab initio) + REvolutionH-tl
  Gene trees visualized on iTOL
```

---

## Repository Structure

```
HBB-Genomics-Pipeline/
|
|-- README.md
|
|-- Task1_GenomicData 
|   |-- Task1_hbb_seq.fasta                      All 7 HBB sequences (FASTA format)
|
|-- Task2_Genome_and_Annotation_Visualization/
|   |-- Task2_hbb.bam                        Coordinate-sorted BAM file
|   |-- Task2_hbb.bai                        BAM index file
|   `-- Task2_igv_session.xml                Saved IGV session (load directly)
|
|-- Task3_Sequence_Alignment_and_Similarity/
|   |-- Task3_ClustalW.clustal               ClustalW 2.1 MSA output (Galaxy)
|   `-- Task3_Muscle.aln-clustalw            MUSCLE 3.8 MSA output (EMBL-EBI)
|
|   [Note: NCBI BLAST pairwise results were web-based; screenshots are in the
|    project report. All 6 runs: Score=2970, Cover=100%, Identity=100%,
|    E-value=0.0, Gaps=0/1608]
|
|-- Task4_Identification_of_Repetitive_Elements/
|   |-- hbb_seq.fasta.out                    Repeat annotations table
|   |-- hbb_seq.fasta.out.gff                GFF3 format (IGV-compatible)
|   |-- hbb_seq.fasta.masked                 Repeat-masked FASTA
|   |-- hbb_seq.fasta.tbl                    RepeatMasker summary statistics
|   |-- hbb_seq.fasta.cat                    Raw alignment data
|   `-- hbb_seq.fasta.out.html               HTML visual summary report
|
|-- Task5_Variant_Calling_and_Functional_Annotation/
|   |-- Task5-FreeBayes_variants.vcf         Raw FreeBayes VCF
|   |-- Task5-VCFfilter.vcf                  Filtered VCF (QUAL>=20, DP>=10)
|   |-- Task5-SnpEff_CSV stats.txt           SnpEff CSV statistics
|   |-- Task5_Map_with_BWA-MEM               Mapped reads in BAM format
|   `-- Task5_SnpEff_HTML_stats/
|       |-- SnpEff_HTML_stats.html           SnpEff HTML report
|       `-- snpeff_stats.genes.txt           Per-gene statistics
|
|-- Task6_Gene_Prediction_and_Orthology_Inference/
|   |-- Augustus_Output.gff                  GFF3 gene predictions (all 7 seqs)
|   |-- Augustus_Predicted_Aminoacidseq.fa   Predicted protein sequences (FASTA)
|   `-- AugustusPredicted_codingseq.fa       Predicted coding sequences (FASTA)
|
|   [Note: REvolutionH-tl was run locally on a four-species bacterial dataset.
|    Output files were generated but are not included here due to size.
|    Full interpretation and screenshots are in the project report.]
|
|-- scripts/
|   |-- 01_fasta_stats.py                    Sequence length and GC content table
|   |-- 02_gc_content_analysis.py            Sliding window GC content analysis
|   |-- 03_repeat_visualizer.py              RepeatMasker results parser and plot
|   |-- 04_vcf_parser.py                     VCF metadata and pipeline summary
|   |-- 05_augustus_parser.py                AUGUSTUS GFF3 parser and protein check
|   `-- 06_pipeline_summary.py               End-to-end pipeline results report
|
`-- figures/                                 Auto-generated figures (see scripts)
    |-- sequence_statistics.png
    |-- gc_content_sliding_window.png
    `-- repeat_distribution.png
```

---

## Tools and Software

| Task | Tool | Version | Platform |
|------|------|---------|---------|
| Data Retrieval | NCBI Nucleotide + BLASTN | — | Web |
| Genome Visualization | IGV Desktop | 2.17+ | Desktop |
| BAM Preparation | BWA-MEM2, SAMtools | — | Galaxy |
| Pairwise Alignment | NCBI BLAST (Megablast) | — | Web |
| Multiple Alignment | ClustalW | 2.1 | Galaxy |
| Multiple Alignment | MUSCLE | 3.8 | EMBL-EBI |
| Repeat Masking | RepeatMasker | 4.2.3 | Local (Ubuntu/WSL2) |
| Repeat Database | Dfam | 3.9 | — |
| Search Engine | RMBlast | 2.14.1+ | Local |
| Read Alignment | BWA-MEM | 0.7.19 | Galaxy |
| Variant Calling | FreeBayes | 1.3.10 | Galaxy |
| Variant Filtering | VCFtools | — | Galaxy |
| Variant Annotation | SnpEff | 5.2 | Galaxy |
| Gene Prediction | AUGUSTUS | — | Web (uni-greifswald.de) |
| Protein Alignment | DIAMOND | 2.1.8 | Local |
| Orthology Inference | REvolutionH-tl | 1.4.6 | Local |
| Tree Visualization | iTOL | — | Web (itol.embl.de) |

---

## Key Results

### Task 1 — Data Retrieval

- 7 HBB sequences retrieved from NCBI GenBank/RefSeq
- GC content: reference 40.11%, population sequences 40.24% (delta = +0.13%)
- Total dataset: 12,552 bp across 7 sequences
- BLASTN confirmed all sequences: 100% identity, 100% coverage, E-value = 0.0

### Task 2 — IGV Genome Visualization

- All 6 reads align to chr11:5,225,464–5,227,071 (GRCh38/hg38)
- NM = 0 (zero mismatches) and MAPQ = 60 (maximum quality) for all samples
- Gene structure confirmed: 3 exons, 2 introns, minus strand orientation
- No sickle cell mutation (HbS, Glu6Val, rs334) detected in any sample
- Uniform coverage depth of 6 across the full 1,608 bp region

### Task 3 — Sequence Alignment

- BLAST (all 6 pairwise runs): Score = 2970, 1608/1608 identities, 0 gaps, E-value = 0.0
- Dot plots: single straight diagonal in all 6 runs, confirming no structural variants
- ClustalW MSA: every position marked `*` (fully conserved) in the 1,608 bp shared region
- MUSCLE MSA: identical conservation — independent algorithmic confirmation
- Conclusion: all 6 population sequences are orthologs of human HBB

### Task 4 — Repeat Masking

- 14 simple repeat elements: (TTTTC)n (7 hits, 252 bp) and (TTTC)n (7 hits, 259 bp)
- Location: positions 639–711 bp on reference (within Intron 2, IVS-2)
- 511 bp masked = 4.07% of total sequence; no transposable elements found
- Masked FASTA (hbb_seq.fasta.masked) generated for downstream use

### Task 5 — Variant Calling

- BWA-MEM, FreeBayes, VCFtools, SnpEff pipeline executed on Galaxy without errors
- FreeBayes VCF: format VCFv4.2, source FreeBayes v1.3.10
- VCFtools filter applied: QUAL >= 20, DP >= 10 (documented in VCF header)
- SnpEff: hg38 loaded (3,247,448,939 bp), 0 warnings, 0 errors
- Zero variants annotated — consistent with wild-type HBB sequences in the dataset

### Task 6 — Gene Prediction and Orthology

- AUGUSTUS predicted 3-exon / 2-intron structure for all 7 sequences (ab initio)
- Gene g1 (reference): span 51–1474 bp, confidence 0.93; Exon 2 confidence 1.0
- Predicted protein (146 aa): exact match to canonical beta-globin (UniProt P68871)
- REvolutionH-tl: all 6 steps completed; gene tree matches Terrabacteria species phylogeny
- HBB confirmed as single-copy vertebrate ortholog in OrthoDB (EOG09150CZZ)

---

## Running the Scripts

### Requirements

```bash
pip install biopython matplotlib pandas seaborn
```

Python 3.8 or later is required. No other installation is needed.

### Usage

Run all scripts from the repository root directory.

```bash
# Sequence length and GC content statistics for all 7 sequences
python scripts/01_fasta_stats.py

# Sliding window GC content across the reference with exon/repeat annotations
python scripts/02_gc_content_analysis.py

# RepeatMasker results summary and position visualization
python scripts/03_repeat_visualizer.py

# VCF file metadata and pipeline step report
python scripts/04_vcf_parser.py

# AUGUSTUS GFF3 gene structure summary and protein comparison
python scripts/05_augustus_parser.py

# Full end-to-end pipeline results report (all tasks)
python scripts/06_pipeline_summary.py
```

Figures are saved automatically to the `figures/` directory.

### Loading the IGV Session

1. Download IGV Desktop: https://igv.org
2. Open IGV and set the genome to **Human (GRCh38/hg38)**
3. File > Load from File > select `Task2_IGV/Task2_hbb.bam`
   (the `.bai` index file must be in the same folder)
4. Navigate to: `chr11:5,225,404-5,227,132`
5. Or: File > Open Session > `Task2_IGV/Task2_igv_session.xml`

---

## Biological Context

The HBB gene encodes the beta-globin chain of adult hemoglobin (HbA), which carries oxygen in red blood cells. A point mutation at codon 6 converting GAG to GTG (Glu to Val, rs334) produces sickle hemoglobin (HbS), which polymerises under low oxygen conditions and deforms red blood cells. This is the molecular basis of Sickle Cell Disease.

The main finding across all six tasks is the **complete conservation** of the HBB coding sequence in all six African population samples. This is consistent with strong purifying selection: because most HBB mutations cause lethal or severely debilitating disease, the functional sequence is maintained across geographically diverse populations.

The (TTTC)n microsatellite found in Intron 2 (IVS-2) by RepeatMasker is a conserved structural feature of this locus, present at the same position in all 7 sequences. Its AT-rich composition is typical of intronic simple sequence repeats and is consistent with the open chromatin environment regulated by the HBB Locus Control Region (LCR).

---

## References

1. NCBI Nucleotide Database. https://www.ncbi.nlm.nih.gov/nucleotide
2. Thorvaldsdottir, H., Robinson, J.T., Mesirov, J.P. (2012). IGV. *Briefings in Bioinformatics*, 14(2), 178–192.
3. Altschul, S.F. et al. (1990). Basic local alignment search tool. *Journal of Molecular Biology*, 215, 403–410.
4. Larkin, M. et al. (2007). ClustalW and ClustalX version 2.0. *Bioinformatics*, 23(21), 2947–2948.
5. Edgar, R.C. (2004). MUSCLE. *Nucleic Acids Research*, 32(5), 1792–1797.
6. Smit, A.F.A., Hubley, R., Green, P. RepeatMasker Open-4.0. http://www.repeatmasker.org
7. Li, H., Durbin, R. (2009). BWA. *Bioinformatics*, 25(14), 1754–1760.
8. Garrison, E., Marth, G. (2012). FreeBayes. *arXiv:1207.3907*
9. Danecek, P. et al. (2011). VCFtools. *Bioinformatics*, 27(15), 2156–2158.
10. Cingolani, P. et al. (2012). SnpEff. *Fly*, 6(2), 80–92.
11. Stanke, M. et al. (2006). AUGUSTUS. *Nucleic Acids Research*, 34, W435–W439.
12. Buchfink, B. et al. (2015). DIAMOND. *Nature Methods*, 12, 59–60.
13. Afgan, E. et al. (2018). Galaxy platform. *Nucleic Acids Research*, 46(W1), W537–W544.
14. Storer, J. et al. (2021). Dfam community resource. *Mobile DNA*, 12(1), 2.

---

## Authors

| Name | Roll No |
|------|---------|
| Namra Basharat | 476203 |
| Ghania Munir | 460673 |
| Hania Fahad | 455132 |
| Nawal Babar | 459508 |
