#!/usr/bin/env python3

import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from pybedtools import BedTool


def load_reference(ref_path):
    record = next(SeqIO.parse(ref_path, "fasta"))
    return str(record.seq), record.id

def load_regions(bed_path):
    return list(BedTool(bed_path).merge())

def get_snps_in_regions(read, ref_seq, bed_regions):
    snps = []

    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=True)

    for qpos, rpos, base in aligned_pairs:
        if qpos is None or rpos is None:
            continue
        # Check if ref position is in any region
        if not any(region.start <= rpos < region.end for region in bed_regions):
            continue

        read_base = base
        ref_base = ref_seq[rpos]

        if read_base != ref_base and read_base != "N" and ref_base != "N":
            snps.append("{}{}{}".format(ref_base, rpos + 1, read_base))

    return snps


def main():
    parser = argparse.ArgumentParser(description = "Summarise SNP combinations in ONT reads within BED regions")
    parser.add_argument("-b", "--bam", required = True, help = "Input BAM file (sorted and indexed)")
    parser.add_argument("-r", "--ref", required = True, help = "Reference FASTA file")
    parser.add_argument("--bed", required = True, help = "BED file with target regions")
    parser.add_argument("-o", "--output", default = "variant_summary.tsv", help = "Output file")
    parser.add_argument("--max-snvs-per-group", type=int, default=None,
                        help="Exclude variant groups with more than this many SNVs (optional)")
    parser.add_argument("--min-frequency", type=int, default=None,
                        help="Exclude variant groups observed fewer than this many times (optional)")

    args = parser.parse_args()

    ref_seq, ref_name = load_reference(args.ref)
    regions = load_regions(args.bed)

    variant_counts = defaultdict(int)

    bamfile = pysam.AlignmentFile(args.bam, "rb")

    skipped_reads = 0

    for read in bamfile.fetch(ref_name):
        if read.is_unmapped or read.query_sequence is None:
            continue

        snps = get_snps_in_regions(read, ref_seq, regions)

        key = ",".join(sorted(snps)) if snps else "ref"
        variant_counts[key] += 1

    bamfile.close()

    with open(args.output, "w") as out:
        for var, count in sorted(variant_counts.items(), key=lambda x: -x[1]):
            if args.max_snvs_per_group is not None and var != "ref":
                snv_count = var.count(",") + 1  # count number of SNVs in group
                if snv_count > args.max_snvs_per_group:
                    continue
            if args.min_frequency is not None and count < args.min_frequency:
                continue
            out.write("{}\t{}\n".format(var, count))



if __name__ == "__main__":
    main()
