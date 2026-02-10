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

# --- Helper for Shell Commands ---
def run_cmd(cmd, description):
    print(f"--- {description} ---")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[!] ERROR: Command failed during: {description}")
        print(f"FAILED COMMAND: {cmd}")
        sys.exit(1)

# --- Helper: Parse BED File (Multi-Region & Length Support) ---
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
                
                # Custom Lengths per Region (Cols 5 & 6)
                r_min = int(parts[4]) if len(parts) >= 6 else default_min
                r_max = int(parts[5]) if len(parts) >= 6 else default_max

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
    ref_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta_path, "fasta"))
    target = None
    
    if chrom in ref_dict:
        target = ref_dict[chrom]
    else:
        # Fuzzy match for 'chr1' vs '1'
        for key in ref_dict:
            if chrom in key or key in chrom:
                target = ref_dict[key]
                chrom = key
                break
    
    if target:
        r_end = min(end, len(target.seq))
        slice_seq = target.seq[start:r_end]
        rec = SeqRecord(slice_seq, id=chrom, description=f"slice_{start}_{end}")
        SeqIO.write(rec, output_path, "fasta")
        return True
    return False

# --- Core Logic: Parse CIGAR for Variants ---
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

# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(description="QONTAS: Hybrid VSEARCH+Minimap Amplicon Pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    parser.add_argument("-b", "--base", required=True, help="Output Basename")
    parser.add_argument("-o", "--outdir", default="Qontas_Out", help="Output Directory")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--region", required=True, help="BED file")
    parser.add_argument("--mincount", type=int, default=10, help="Minimum sequence count (denoising threshold)")
    parser.add_argument("--ignore_edges", type=int, default=3, help="Ignore variants at slice ends")
    parser.add_argument("--default_minlen", type=int, default=50)
    parser.add_argument("--default_maxlen", type=int, default=10000)

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    regions = parse_bed_regions(args.region, args.default_minlen, args.default_maxlen)
    
    # 1. Global Initial Alignment
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    if not os.path.exists(full_bam):
        run_cmd(f"minimap2 -a -x map-ont -t {args.threads} {args.ref} {args.input} | "
                f"samtools view -u -F 2308 | samtools sort -o {full_bam}", "Global Alignment")
        run_cmd(f"samtools index {full_bam}", "Indexing Full BAM")

    bam = pysam.AlignmentFile(full_bam, "rb")

    for reg in regions:
        prefix = f"{args.base}_{reg['name']}"
        print(f"\n>>> Processing Region: {reg['name']}")
        
        # A. Extract Reads for Region
        reg_fa = os.path.join(args.outdir, f"{prefix}_temp.fa")
        count = 0
        with open(reg_fa, "w") as f_out:
            try:
                for r in bam.fetch(reg['chrom'], reg['start'], reg['end']):
                    # Coverage Check: Read must span the whole region
                    if r.reference_start <= reg['start'] and r.reference_end >= reg['end']:
                        if reg['min_len'] <= r.query_length <= reg['max_len']:
                            f_out.write(f">{r.query_name}\n{r.query_sequence}\n")
                            count += 1
            except ValueError:
                print(f"    [!] Chromosome {reg['chrom']} not found in BAM.")
                continue
        
        if count < args.mincount:
            print(f"    [!] Insufficient reads ({count}). Skipping.")
            if os.path.exists(reg_fa): os.remove(reg_fa)
            continue

        # B. Denoise with VSEARCH
        reg_derep = os.path.join(args.outdir, f"{prefix}_derep.fa")
        reg_uc = os.path.join(args.outdir, f"{prefix}_derep.txt")
        run_cmd(f"vsearch --derep_fulllength {reg_fa} --output {reg_derep} --uc {reg_uc} "
                f"--minuniquesize {args.mincount} --sizeout", "VSEARCH Denoising")

        # C. Map Unique Sequences (Dynamic Preset Selection)
        ref_slice = os.path.join(args.outdir, f"{prefix}_ref_slice.fa")
        if not create_ref_slice_fasta(args.ref, reg['chrom'], reg['start'], reg['end'], ref_slice): 
            continue
        
        region_len = reg['end'] - reg['start']
        preset = "sr" if region_len < 300 else "map-ont"
        unique_bam = os.path.join(args.outdir, f"{prefix}_unique.bam")
        
        run_cmd(f"minimap2 -a --eqx -x {preset} -t {args.threads} {ref_slice} {reg_derep} | "
                f"samtools sort -o {unique_bam}", f"Final Mapping (Preset: {preset})")
        run_cmd(f"samtools index {unique_bam}", "Indexing Unique BAM")

        # D. Parse CIGARs
        ref_slice_str = str(next(SeqIO.parse(ref_slice, "fasta")).seq)
        var_map = {}
        with pysam.AlignmentFile(unique_bam, "rb") as ubam:
            mapped_count = 0
            for r in ubam.fetch():
                if not (r.is_secondary or r.is_supplementary or r.is_unmapped):
                    mapped_count += 1
                    # Harmonize ID: strip size metadata for matching
                    clean_id = r.query_name.split(';')[0]
                    var_map[clean_id] = parse_cigar_for_variants(r, ref_slice_str, args.ignore_edges)
        
        # E. Table Generation
        try:
            df_uc = pd.read_csv(reg_uc, sep="\t", header=None, dtype=str)
            # 'C' rows contain cluster summaries
            clusters = df_uc[df_uc[0] == "C"][[8, 2]].copy()
            clusters.columns = ["Cluster", "Count"]
            clusters["Cluster"] = clusters["Cluster"].str.split(';').str[0]
            clusters["Count"] = pd.to_numeric(clusters["Count"])
            
            # Match variants from the map
            clusters["Variant"] = clusters["Cluster"].map(var_map).fillna("Unaligned")
            
            # Summarize by Variant Signature
            final = clusters.groupby("Variant")["Count"].sum().reset_index()
            final["RelAbund"] = (final["Count"] / final["Count"].sum()) * 100
            final = final.sort_values("Count", ascending=False)
            
            out_tsv = os.path.join(args.outdir, f"{prefix}_variants.tsv")
            final.to_csv(out_tsv, sep="\t", index=False)
            print(f"    [Debug] Mapped {mapped_count} clusters. Results saved to {out_tsv}")
        except Exception as e:
            print(f"    [!] Error generating table: {e}")

        # Final Cleanup for Region
        for f in [reg_fa, reg_derep, reg_uc, ref_slice, unique_bam, unique_bam+".bai"]:
            if os.path.exists(f): os.remove(f)

    bam.close()
    print(f"\nSUCCESS. Final files are in {args.outdir}")

if __name__ == "__main__":
    main()