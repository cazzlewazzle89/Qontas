#!/usr/bin/env python3

import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def write_forward_reads(input_bam, output_fasta):
    bam = pysam.AlignmentFile(input_bam, "rb")

    forward_reads = []
    for read in bam.fetch():
        header = read.query_name
        seq = read.query_sequence
        forward_reads.append((header, seq))

    bam.close()

    seq_records = []
    for header, sequence in forward_reads:
        seq = Seq(sequence)
        record = SeqRecord(seq, id = header)
        seq_records.append(record)

    SeqIO.write(seq_records, output_fasta, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write forward-oriented reads from a BAM file to a FASTA file")
    parser.add_argument("--input_bam", dest='input_bam', help="Input BAM file path")
    parser.add_argument("--output_fasta", dest='output_fasta', help="Output FASTA file path")
    args = parser.parse_args()

    write_forward_reads(args.input_bam, args.output_fasta)