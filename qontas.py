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

    args = parser.parse_args()

    ref_seq, ref_name = load_reference(args.ref)
    regions = load_regions(args.bed)

    variant_counts = defaultdict(int)

    bamfile = pysam.AlignmentFile(args.bam, "rb")

    for read in bamfile.fetch(ref_name):
        if read.is_unmapped:
            continue

        # Get aligned reference positions
        ref_pos = read.get_reference_positions(full_length=True)
        seq = read.query_sequence

        if read.is_reverse:
            seq = str(Seq(seq).reverse_complement())
            ref_pos = ref_pos[::-1]

        read.query_sequence = seq  

        snps = get_snps_in_regions(read, ref_seq, regions)

        key = ",".join(sorted(snps)) if snps else "ref"
        variant_counts[key] += 1

    bamfile.close()

    with open(args.output, "w") as out:
        for var, count in sorted(variant_counts.items(), key=lambda x: -x[1]):
            out.write("{}\t{}\n".format(var, count))


if __name__ == "__main__":
    main()
