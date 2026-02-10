#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd, description):
    print(f"--- {description} ---")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[!] ERROR: Command failed during: {description}")
        print(f"FAILED COMMAND: {cmd}")
        sys.exit(1)

def parse_bed_regions(bed_file, default_min, default_max):
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): continue
            parts = line.split('\t')
            if len(parts) < 3: continue
            chrom, start, end = parts[0].strip(), int(parts[1]), int(parts[2])
            name = parts[3].strip() if len(parts) >= 4 else f"reg_{parts[0]}_{start}"
            r_min = int(parts[4]) if len(parts) >= 6 else default_min
            r_max = int(parts[5]) if len(parts) >= 6 else default_max
            regions.append({'chrom': chrom, 'name': name, 'start': start, 'end': end, 'min_len': r_min, 'max_len': r_max})
    return regions

def parse_cigar_for_variants(read, ref_seq_str, ignore_edges=3):
    if not read.cigartuples: return "Unaligned"
    ref_pos, query_pos = 0, 0
    variants = []
    ref_len = len(ref_seq_str)
    
    for op, length in read.cigartuples:
        # 7=EQ (=), 8=X (mismatch), 1=I, 2=D, 4=S, 0=M
        if op == 7: 
            ref_pos += length
            query_pos += length
        elif op == 8:
            for k in range(length):
                if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                    r_base = ref_seq_str[ref_pos]
                    q_base = read.query_sequence[query_pos + k]
                    variants.append(f"{r_base}{ref_pos + 1}{q_base}")
                ref_pos += 1
            query_pos += length
        elif op == 1:
            if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                ins_seq = read.query_sequence[query_pos : query_pos + length]
                variants.append(f"ins{ref_pos}{ins_seq}")
            query_pos += length
        elif op == 2:
            for k in range(length):
                if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                    variants.append(f"del{ref_pos + 1}")
                ref_pos += 1
        elif op in [0, 4]:
            ref_pos += length if op == 0 else 0
            query_pos += length
    return ",".join(variants) if variants else "Ref"

def main():
    parser = argparse.ArgumentParser(description="QONTAS: Alignment-First Binning & Quantification")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-b", "--base", required=True)
    parser.add_argument("-o", "--outdir", default="Qontas_Out")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--region", required=True)
    parser.add_argument("--mincount", type=int, default=10)
    parser.add_argument("--ignore_edges", type=int, default=3)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    regions = parse_bed_regions(args.region, 50, 10000)
    
    # 1. Global Alignment
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    if not os.path.exists(full_bam):
        run_cmd(f"minimap2 -a -x map-ont -t {args.threads} {args.ref} {args.input} | samtools view -u -F 2308 | samtools sort -o {full_bam}", "Global Alignment")
        run_cmd(f"samtools index {full_bam}", "Indexing Full BAM")

    bam = pysam.AlignmentFile(full_bam, "rb")

    for reg in regions:
        prefix = f"{args.base}_{reg['name']}"
        print(f"\n>>> Processing: {reg['name']}")
        
        # 2, 3, 4. Bin, Filter, and Hard Clip
        reg_fa = os.path.join(args.outdir, f"{prefix}_clipped.fa")
        count = 0
        with open(reg_fa, "w") as f_out:
            try:
                for r in bam.fetch(reg['chrom'], reg['start'], reg['end']):
                    # Check Coverage
                    if r.reference_start <= reg['start'] and r.reference_end >= reg['end']:
                        # Check Length
                        if reg['min_len'] <= r.query_length <= reg['max_len']:
                            
                            # Hard Clip to BED region
                            q_start, q_end = None, None
                            aligned_pairs = r.get_aligned_pairs(matches_only=True)
                            for q_p, r_p in aligned_pairs:
                                if q_start is None and r_p >= reg['start']: q_start = q_p
                                if r_p <= reg['end']: q_end = r_p
                            
                            if q_start is not None and q_end is not None:
                                s, e = min(q_start, q_end), max(q_start, q_end)
                                # Extract and Reorient
                                seq_obj = Seq(r.query_sequence[s:e])
                                if r.is_reverse:
                                    seq_obj = seq_obj.reverse_complement()
                                
                                f_out.write(f">{r.query_name}\n{str(seq_obj)}\n")
                                count += 1
            except ValueError: continue

        if count < args.mincount:
            print(f"    [!] Insufficient reads ({count}). Skipping.")
            continue

        # 5. Count Unique Sequences (Denoise)
        reg_derep = os.path.join(args.outdir, f"{prefix}_derep.fa")
        reg_uc = os.path.join(args.outdir, f"{prefix}_derep.txt")
        run_cmd(f"vsearch --derep_fulllength {reg_fa} --output {reg_derep} --uc {reg_uc} --minuniquesize {args.mincount} --sizeout", "Denoising")

        # 6. Name Sequences via Minimap2 CIGAR
        ref_slice_fa = os.path.join(args.outdir, f"{prefix}_ref_slice.fa")
        ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))
        target_seq = ref_dict[reg['chrom']]
        slice_rec = SeqRecord(target_seq.seq[reg['start']:reg['end']], id=reg['chrom'], description="slice")
        SeqIO.write(slice_rec, ref_slice_fa, "fasta")
        
        region_len = reg['end'] - reg['start']
        preset = "sr" if region_len < 300 else "map-ont"
        unique_bam = os.path.join(args.outdir, f"{prefix}_unique.bam")
        run_cmd(f"minimap2 -a --eqx -x {preset} -t {args.threads} {ref_slice_fa} {reg_derep} | samtools sort -o {unique_bam}", "Naming Alignment")
        run_cmd(f"samtools index {unique_bam}", "Indexing Unique BAM")

        # Assemble Results
        ref_slice_str = str(slice_rec.seq)
        var_map = {}
        with pysam.AlignmentFile(unique_bam, "rb") as ubam:
            for r in ubam.fetch():
                if not (r.is_secondary or r.is_supplementary or r.is_unmapped):
                    clean_id = r.query_name.split(';')[0]
                    var_map[clean_id] = parse_cigar_for_variants(r, ref_slice_str, args.ignore_edges)

        df = pd.read_csv(reg_uc, sep="\t", header=None, dtype=str)
        clusters = df[df[0] == "C"][[8, 2]].copy()
        clusters.columns = ["Cluster", "Count"]
        clusters["Cluster_Clean"] = clusters["Cluster"].str.split(';').str[0]
        clusters["Variant"] = clusters["Cluster_Clean"].map(var_map).fillna("Unaligned")

        final = clusters.groupby("Variant")["Count"].sum().reset_index().sort_values("Count", ascending=False)
        final["RelAbund"] = (final["Count"] / final["Count"].sum()) * 100
        final.to_csv(os.path.join(args.outdir, f"{prefix}_variants.tsv"), sep="\t", index=False)
        
        for f in [reg_fa, reg_derep, reg_uc, ref_slice_fa, unique_bam, unique_bam+".bai"]:
            if os.path.exists(f): os.remove(f)

    print(f"\nSUCCESS. Check {args.outdir}")

if __name__ == "__main__":
    main()
    