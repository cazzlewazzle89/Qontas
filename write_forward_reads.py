#!/usr/bin/env python3

import argparse
import pysam
from Bio import SeqIO

def write_forward_reads(bam_file, output_fasta):
    bam = pysam.AlignmentFile(bam_file, "rb")

    forward_reads = []
    for read in bam.fetch():
        header = read.query_name
        seq = read.query_sequence
        forward_reads.append((header, seq))

    bam.close()

    with open(output_fasta, "w") as fasta_file:
        for header, seq in forward_reads:
            fasta_file.write(f">{header}\n{seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write forward-oriented reads from a BAM file to a FASTA file")
    parser.add_argument("--bam_file", dest='bam_file', help="Input BAM file path")
    parser.add_argument("--output_fasta", dest='output_fasta', help="Output FASTA file path")
    args = parser.parse_args()

    write_forward_reads(args.bam_file, args.output_fasta)