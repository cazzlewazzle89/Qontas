#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import gzip
import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- Helper for Shell Commands ---
def run_cmd(cmd, description):
    print(f"--- {description} ---")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[!] ERROR: Command failed during: {description}")
        print(f"FAILED COMMAND: {cmd}")
        sys.exit(1)

# --- 1. Extract & Trim Reads (Hard Clipping + Reorientation) ---
def extract_and_trim_reads(input_bam, output_fasta, region_bed=None):
    print("--- Extracting and Trimming Reads ---")
    
    ref_start, ref_end = None, None
    if region_bed:
        try:
            with open(region_bed, 'r') as f:
                # Assuming simple BED: Chrom Start End
                parts = f.readline().strip().split('\t')
                ref_start = int(parts[1])
                ref_end = int(parts[2])
            print(f"    Trimming to Target Region: {ref_start}-{ref_end}")
        except Exception as e:
            print(f"Error parsing BED file: {e}")
            sys.exit(1)

    bam = pysam.AlignmentFile(input_bam, "rb")
    seq_records = []
    processed = 0
    skipped = 0

    for read in bam.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue

        # REGION FILTERING
        if ref_start is not None and ref_end is not None:
            # Skip if read doesn't FULLY cover region
            if read.reference_start > ref_start or read.reference_end < ref_end:
                skipped += 1
                continue

            q_start = None
            q_end = None
            
            # Robust CIGAR mapping
            aligned_pairs = read.get_aligned_pairs(matches_only=True)
            for q_pos, r_pos in aligned_pairs:
                if q_start is None and r_pos >= ref_start:
                    q_start = q_pos
                if r_pos <= ref_end:
                    q_end = q_pos

            if q_start is not None and q_end is not None:
                start_idx = min(q_start, q_end)
                end_idx = max(q_start, q_end)
                
                # Slicing
                seq_str = read.query_sequence[start_idx:end_idx]
                seq_obj = Seq(seq_str)
                
                # Reorienting
                if read.is_reverse:
                    seq_obj = seq_obj.reverse_complement()

                record = SeqRecord(seq_obj, id=read.query_name, description="")
                seq_records.append(record)
                processed += 1
            else:
                skipped += 1
        
        else:
            # Global extraction (no trimming)
            seq_obj = Seq(read.query_sequence)
            if read.is_reverse:
                seq_obj = seq_obj.reverse_complement()
            record = SeqRecord(seq_obj, id=read.query_name, description="")
            seq_records.append(record)
            processed += 1

    bam.close()
    SeqIO.write(seq_records, output_fasta, "fasta")
    print(f"    Extracted {processed} reads. Skipped {skipped}.")

# --- 2. Make Feature Table (Corrected Logic) ---
def make_feature_table(input_uc, derep_fasta, output_prefix, min_count, min_abund):
    print("--- Generating Feature Table ---")
    
    # Load vsearch UC file
    try:
        df = pd.read_csv(input_uc, sep="\t", header=None, on_bad_lines='skip')
    except Exception:
        df = pd.read_csv(input_uc, sep="\t", header=None, error_bad_lines=False)

    # Filter for Cluster ('C') rows only
    # Col 8 = Cluster ID, Col 2 = Cluster Size
    clusters = df[df[0] == "C"][[8, 2]].copy()
    clusters.columns = ["Cluster", "Count"]

    # Filter by Count
    clusters = clusters[clusters["Count"] >= min_count]
    
    # Calculate Relative Abundance
    total_count = clusters["Count"].sum()
    if total_count > 0:
        clusters["RelAbund"] = 100 * clusters["Count"] / total_count
        # Filter by Abundance
        clusters = clusters[clusters["RelAbund"] >= min_abund]
    else:
        clusters["RelAbund"] = 0.0

    # Output 1: TSV Table
    clusters_output = f"{output_prefix}_clusters.txt"
    clusters.to_csv(clusters_output, sep="\t", index=False, header=False)
    print(f"    Table written to {clusters_output}")

    # Output 2: Readnames list
    readnames_output = f"{output_prefix}_readnames.txt"
    with open(readnames_output, 'w') as f:
        for text in clusters['Cluster'].tolist():
            f.write(str(text) + '\n')

    # Output 3: Clean FASTA
    # We load the full set of unique sequences (derep_fasta) and keep only those in our table
    valid_ids = set(clusters['Cluster'].astype(str))
    
    final_seqs_file = f"{output_prefix}_unique_sequences.fasta"
    count_seqs = 0
    with open(final_seqs_file, "w") as out_f:
        for record in SeqIO.parse(derep_fasta, "fasta"):
            if record.id in valid_ids:
                SeqIO.write(record, out_f, "fasta")
                count_seqs += 1
                
    print(f"    Filtered Unique Sequences written to {final_seqs_file} ({count_seqs} seqs).")

# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(description="QONTAS: Microbiome Amplicon Pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    parser.add_argument("-b", "--base", required=True, help="Output Basename")
    parser.add_argument("-o", "--outdir", required=True, help="Output Directory")
    parser.add_argument("--minlen", type=int, default=600)
    parser.add_argument("--maxlen", type=int, default=650)
    parser.add_argument("--mincount", type=int, default=10)
    parser.add_argument("--minabund", type=float, default=1.0)
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--region", help="Optional BED file for trimming sequences")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Filter FASTQ (Using SeqKit for multithreading speed)
    filt_fastq = os.path.join(args.outdir, f"{args.base}_filtered.fastq.gz")
    
    # Check if seqkit is installed/available
    seqkit_cmd = (f"seqkit seq -j {args.threads} "
                  f"-m {args.minlen} -M {args.maxlen} "
                  f"{args.input} -o {filt_fastq}")
    run_cmd(seqkit_cmd, "Filtering FASTQ with SeqKit")

    # 2. Map and Sort
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    map_cmd = (f"minimap2 -a -x map-ont -t {args.threads} {args.ref} {filt_fastq} | "
               f"samtools view -u -F 2308 | "
               f"samtools sort -@ {args.threads} -o {full_bam}")
    run_cmd(map_cmd, "Mapping with Minimap2")
    run_cmd(f"samtools index {full_bam}", "Indexing BAM")

    # 3. Extract & Trim
    final_fasta = os.path.join(args.outdir, f"{args.base}.fasta")
    extract_and_trim_reads(full_bam, final_fasta, args.region)

    # 4. Dereplicate
    derep_fasta = os.path.join(args.outdir, f"{args.base}_derep.fasta")
    derep_txt = os.path.join(args.outdir, f"{args.base}_derep.txt")
    run_cmd(f"vsearch --derep_fulllength {final_fasta} --output {derep_fasta} --threads {args.threads} --uc {derep_txt}",
            "Dereplicating with VSEARCH")

    # 5. Feature Table & Unique Seqs
    out_prefix = os.path.join(args.outdir, args.base)
    make_feature_table(derep_txt, derep_fasta, out_prefix, args.mincount, args.minabund)

    # Cleanup
    print("\n--- Cleaning up ---")
    if os.path.exists(full_bam): os.remove(full_bam)
    if os.path.exists(full_bam + ".bai"): os.remove(full_bam + ".bai")
    # if os.path.exists(filt_fastq): os.remove(filt_fastq) # Keep if you want to inspect filtered reads

    print(f"\nSUCCESS: Pipeline complete. Results in {args.outdir}")

if __name__ == "__main__":
    main()