#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
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

# --- Helper: Parse BED File ---
def parse_bed_regions(bed_file, default_min, default_max):
    regions = []
    try:
        with open(bed_file, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split('\t')
                if len(parts) < 3: continue
                
                chrom = parts[0].strip()
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3].strip() if len(parts) >= 4 else f"region_{i+1:02d}"
                
                if len(parts) >= 6:
                    r_min, r_max = int(parts[4]), int(parts[5])
                else:
                    r_min, r_max = default_min, default_max

                regions.append({
                    'chrom': chrom, 'name': name, 'start': start, 'end': end,
                    'min_len': r_min, 'max_len': r_max
                })
    except Exception as e:
        print(f"Error parsing BED file: {e}")
        sys.exit(1)
    return regions

# --- Helper: Create Reference Slice File ---
def create_ref_slice_fasta(ref_fasta_path, chrom, start, end, output_path):
    records = list(SeqIO.parse(ref_fasta_path, "fasta"))
    target = None
    for r in records:
        if r.id.strip() == chrom:
            target = r
            break
    if target is None and len(records) == 1:
        target = records[0]

    if not target:
        return False
        
    r_end = min(end, len(target.seq))
    slice_seq = target.seq[start:r_end]
    rec = SeqRecord(slice_seq, id=chrom, description="slice")
    SeqIO.write(rec, output_path, "fasta")
    return True

# --- Core Logic: Parse CIGAR for Variants ---
def parse_cigar_for_variants(cigar_tuples, query_seq, ref_seq_str, ignore_edges=3):
    ref_pos = 0 
    query_pos = 0
    variants = []
    ref_len = len(ref_seq_str)
    
    for op, length in cigar_tuples:
        if op == 7: # EQ
            ref_pos += length
            query_pos += length
        elif op == 8: # X
            for k in range(length):
                if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                    r_base = ref_seq_str[ref_pos]
                    q_base = query_seq[query_pos + k]
                    variants.append(f"{r_base}{ref_pos + 1}{q_base}")
                ref_pos += 1
            query_pos += length
        elif op == 1: # I
            if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                inserted_bases = query_seq[query_pos : query_pos + length]
                variants.append(f"ins{ref_pos}{inserted_bases}")
            query_pos += length
        elif op == 2: # D
            for k in range(length):
                if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                    variants.append(f"del{ref_pos + 1}")
                ref_pos += 1
        elif op == 4: # S
            query_pos += length
        elif op == 0: # M
            ref_pos += length
            query_pos += length

    return ",".join(variants) if variants else "Ref"

# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(description="QONTAS: HQuantification of ONT Amplicon Sequences")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    parser.add_argument("-b", "--base", required=True, help="Output Basename")
    parser.add_argument("-o", "--outdir", default="Qontas_Out", help="Output Directory")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--region", required=True, help="BED file")
    parser.add_argument("--mincount", type=int, default=10)
    parser.add_argument("--ignore_edges", type=int, default=3)
    parser.add_argument("--default_minlen", type=int, default=600)
    parser.add_argument("--default_maxlen", type=int, default=650)

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    regions = parse_bed_regions(args.region, args.default_minlen, args.default_maxlen)
    
    # 1. Map Everything (Global)
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    if not os.path.exists(full_bam):
        map_cmd = (f"seqkit seq -m 20 {args.input} | " 
                   f"minimap2 -a -x map-ont -t {args.threads} {args.ref} - | "
                   f"samtools view -u -F 2308 | "
                   f"samtools sort -@ {args.threads} -o {full_bam}")
        run_cmd(map_cmd, "Initial Mapping")
        run_cmd(f"samtools index {full_bam}", "Indexing Initial BAM")

    import pysam
    bam = pysam.AlignmentFile(full_bam, "rb")

    for reg in regions:
        chrom, name, start, end = reg['chrom'], reg['name'], reg['start'], reg['end']
        min_l, max_l = reg['min_len'], reg['max_len']
        
        print(f"--- Processing Region: {name} ---")
        prefix = f"{args.base}_{name}"
        
        # A. Extract
        reg_fasta = os.path.join(args.outdir, f"{prefix}_temp.fasta")
        processed_count = 0
        with open(reg_fasta, 'w') as f_out:
            try:
                for read in bam.fetch(contig=chrom, start=start, stop=end):
                    if read.is_unmapped: continue
                    if read.reference_start > start or read.reference_end < end: continue
                    if read.query_length < min_l or read.query_length > max_l: continue
                    f_out.write(f">{read.query_name}\n{read.query_sequence}\n")
                    processed_count += 1
            except ValueError: continue
        
        if processed_count == 0:
            if os.path.exists(reg_fasta): os.remove(reg_fasta)
            continue

        # B. Denoise (VSEARCH)
        reg_derep = os.path.join(args.outdir, f"{prefix}_derep.fasta")
        reg_uc = os.path.join(args.outdir, f"{prefix}_derep.txt")
        run_cmd(f"vsearch --derep_fulllength {reg_fasta} --output {reg_derep} --uc {reg_uc} --minuniquesize {args.mincount} --sizeout", "Denoising with VSEARCH")

        # C. Map Unique & Index
        ref_slice_fa = os.path.join(args.outdir, f"{prefix}_ref_slice.fasta")
        if not create_ref_slice_fasta(args.ref, chrom, start, end, ref_slice_fa): continue
        
        unique_bam = os.path.join(args.outdir, f"{prefix}_unique.bam")
        # Ensure --eqx for variant calling
        map_u_cmd = (f"minimap2 -a --eqx -x map-ont -t {args.threads} {ref_slice_fa} {reg_derep} | "
                     f"samtools view -b -o {unique_bam}")
        run_cmd(map_u_cmd, "Mapping Unique Sequences")
        run_cmd(f"samtools index {unique_bam}", "Indexing Unique BAM") # FIXED: Added index step
        
        # D. Parse CIGARs
        ref_slice_str = str(next(SeqIO.parse(ref_slice_fa, "fasta")).seq)
        ubam = pysam.AlignmentFile(unique_bam, "rb")
        variant_map = {}
        
        for read in ubam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            var_str = parse_cigar_for_variants(read.cigartuples, read.query_sequence, ref_slice_str, args.ignore_edges)
            variant_map[read.query_name] = var_str
        ubam.close()
        
        # E. Table
        try:
            df = pd.read_csv(reg_uc, sep="\t", header=None, dtype=str)
            clusters = df[df[0] == "C"][[8, 2]].copy()
            clusters.columns = ["Cluster", "Count"]
            clusters["Count"] = pd.to_numeric(clusters["Count"])
            clusters["Variant"] = clusters["Cluster"].map(variant_map).fillna("Filtered/Unaligned")
            
            final_df = clusters.groupby("Variant")["Count"].sum().reset_index()
            final_df["RelAbund"] = (final_df["Count"] / final_df["Count"].sum()) * 100
            final_df = final_df.sort_values(by="Count", ascending=False)
            
            out_tsv = os.path.join(args.outdir, f"{prefix}_variants.tsv")
            final_df.to_csv(out_tsv, sep="\t", index=False)
        except Exception as e:
            print(f"Error processing table: {e}")

        # Cleanup
        for f in [ref_slice_fa, unique_bam, unique_bam+".bai", reg_derep, reg_uc, reg_fasta]:
            if os.path.exists(f): os.remove(f)

    bam.close()
    print(f"\nSUCCESS. Results in {args.outdir}")

if __name__ == "__main__":
    main()
