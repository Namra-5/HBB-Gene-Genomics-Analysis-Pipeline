#!/usr/bin/env python3
"""
Script 04: VCF File Parser and Summary
HBB Genomics Project — Task 5 Bonus Analysis

Parses the FreeBayes raw VCF and VCFtools-filtered VCF files,
extracts metadata and variant records, and produces a structured report.

Usage:
    python 04_vcf_parser.py
    python 04_vcf_parser.py --raw path/to/raw.vcf --filtered path/to/filtered.vcf
"""

import os
import sys

DEFAULT_RAW      = os.path.join(os.path.dirname(__file__), "..", "Task5_Variant_Calling_and_Functional_Annotation",
                                "Task5-FreeBayes_variants.vcf")
DEFAULT_FILTERED = os.path.join(os.path.dirname(__file__), "..", "Task5_Variant_Calling_and_Functional_Annotation",
                                "Task5-VCFfilter.vcf")


def parse_vcf(filepath):
    """
    Parse a VCF file. Returns (meta, header_fields, variants).
    meta         : list of ##... lines
    header_fields: list of column names from the #CHROM line
    variants     : list of dicts (one per data line)
    """
    meta          = []
    header_fields = []
    variants      = []

    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#CHROM"):
                header_fields = line[1:].split("\t")
            else:
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 8:
                    v = {
                        "CHROM":  parts[0],
                        "POS":    parts[1],
                        "ID":     parts[2],
                        "REF":    parts[3],
                        "ALT":    parts[4],
                        "QUAL":   parts[5],
                        "FILTER": parts[6],
                        "INFO":   parts[7],
                    }
                    if len(parts) >= 10:
                        v["FORMAT"] = parts[8]
                        v["SAMPLE"] = parts[9]
                    variants.append(v)

    return meta, header_fields, variants


def extract_meta_value(meta, key):
    """Extract a value from a ##key=value metadata line."""
    for line in meta:
        if line.startswith(f"##{key}="):
            return line.split("=", 1)[1]
    return "not found"


def info_field(info_str, key):
    """Extract a specific key's value from an INFO field string."""
    for part in info_str.split(";"):
        if part.startswith(f"{key}="):
            return part.split("=", 1)[1]
    return "."


def print_report(raw_path, filtered_path):
    print(f"\n{'='*72}")
    print("  VCF File Analysis — Task 5: Variant Calling and Annotation")
    print(f"{'='*72}")

    for label, path in [("RAW VCF (FreeBayes output)", raw_path),
                         ("FILTERED VCF (VCFtools output)", filtered_path)]:
        if not os.path.exists(path):
            print(f"\n  [{label}] File not found: {path}")
            continue

        meta, header, variants = parse_vcf(path)
        filesize = os.path.getsize(path)

        print(f"\n  [{label}]")
        print(f"  File: {os.path.basename(path)}")
        print(f"  Size: {filesize:,} bytes")
        print(f"  Metadata lines (##): {len(meta)}")
        print(f"  Variant records    : {len(variants)}")

        # Key metadata
        print(f"\n  Key Metadata:")
        print(f"    VCF format version : {extract_meta_value(meta, 'fileformat')}")
        print(f"    Date               : {extract_meta_value(meta, 'fileDate')}")
        print(f"    Source tool        : {extract_meta_value(meta, 'source')}")
        print(f"    Reference          : {extract_meta_value(meta, 'reference')}")

        # Contig lines
        contigs = [m for m in meta if m.startswith("##contig")]
        print(f"\n  Contigs defined ({len(contigs)}):")
        for c in contigs:
            # parse ID= and length= from the contig line
            c_inner = c.replace("##contig=<", "").rstrip(">")
            fields  = {f.split("=")[0]: f.split("=")[1] for f in c_inner.split(",") if "=" in f}
            print(f"    ID={fields.get('ID','?'):<40}  length={fields.get('length','?')} bp")

        # Filter lines (present in filtered VCF)
        filter_lines = [m for m in meta if m.startswith("##filter")]
        if filter_lines:
            print(f"\n  Applied Filters:")
            for fl in filter_lines:
                print(f"    {fl}")

        # Variants
        if not variants:
            print(f"\n  Variant Records: 0 data lines found in this VCF.")
            print( "    This is expected for the SnpEff annotation step when no variants")
            print( "    passed all filters within the coordinate system of the full hg38 genome.")
        else:
            print(f"\n  Variant Records ({len(variants)}):")
            print(f"  {'CHROM':<42} {'POS':>6} {'REF':<6} {'ALT':<10} {'QUAL':>7} {'TYPE'}")
            print(f"  {'-'*42} {'-'*6} {'-'*6} {'-'*10} {'-'*7} {'-'*10}")
            for v in variants:
                vtype = info_field(v["INFO"], "TYPE")
                print(f"  {v['CHROM']:<42} {v['POS']:>6} {v['REF']:<6} {v['ALT']:<10} "
                      f"{v['QUAL']:>7} {vtype}")

    # Pipeline summary
    print(f"\n  Pipeline Step Summary:")
    print(f"  {'Step':<30} {'Tool':<18} {'Output':<35} {'Result'}")
    print(f"  {'-'*30} {'-'*18} {'-'*35} {'-'*30}")
    steps = [
        ("1. Read Alignment",       "BWA-MEM 0.7.19", "BAM file (1.6 KB)",           "Alignment complete, hg38"),
        ("2. Variant Calling",      "FreeBayes 1.3.10","Raw VCF (7.8 KB, 8 records)", "QUAL+depth recorded"),
        ("3. Variant Filtering",    "VCFtools",        "Filtered VCF (QUAL>=20, DP>=10)","Filter header added"),
        ("4. Functional Annotation","SnpEff 5.2",      "Annotated VCF + HTML + CSV",  "hg38 DB loaded, 0 errors"),
    ]
    for step, tool, output, result in steps:
        print(f"  {step:<30} {tool:<18} {output:<35} {result}")

    print(f"\n  Note on Zero Variants in SnpEff Output:")
    print(f"    The raw FreeBayes VCF records variants against a local HBB reference")
    print(f"    (1,608 bp contig). SnpEff annotates against the full hg38 genome. When")
    print(f"    variant positions in the local contig coordinate system do not map to a")
    print(f"    recognised chromosome in hg38, SnpEff reports 0 annotated variants.")
    print(f"    This is a coordinate system mismatch, not a pipeline failure.")
    print(f"    The correct approach for production analysis would be to call variants")
    print(f"    with chr11 as the reference, so positions map directly to hg38.")
    print(f"\n{'='*72}\n")


def main():
    raw_path      = DEFAULT_RAW
    filtered_path = DEFAULT_FILTERED

    if "--raw" in sys.argv:
        raw_path = sys.argv[sys.argv.index("--raw") + 1]
    if "--filtered" in sys.argv:
        filtered_path = sys.argv[sys.argv.index("--filtered") + 1]

    # Fallback paths (run from repo root)
    if not os.path.exists(raw_path):
        raw_path = os.path.join("Task5_Variant_Calling_and_Functional_Annotation", "Task5-FreeBayes_variants.vcf")
    if not os.path.exists(filtered_path):
        filtered_path = os.path.join("Task5_Variant_Calling_and_Functional_Annotation", "Task5-VCFfilter.vcf")

    print_report(raw_path, filtered_path)


if __name__ == "__main__":
    main()
