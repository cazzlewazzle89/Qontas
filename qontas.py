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
    # Load Ref (optimized for single lookups)
    records = list(SeqIO.parse(ref_fasta_path, "fasta"))
    target = None
    
    # Try exact match
    for r in records:
        if r.id.strip() == chrom:
            target = r
            break
    # Fallback
    if target is None and len(records) == 1:
        target = records[0]

    if not target:
        print(f"    [!] Error: Chrom '{chrom}' not found in Ref.")
        return False
        
    # Extract Slice (0-based)
    # Ensure we don't go out of bounds
    r_end = min(end, len(target.seq))
    slice_seq = target.seq[start:r_end]
    
    # Write to temp file
    rec = SeqRecord(slice_seq, id=chrom, description="slice")
    SeqIO.write(rec, output_path, "fasta")
    return True

# --- Core Logic: Parse CIGAR for Variants ---
def parse_cigar_for_variants(cigar_tuples, query_seq, ref_seq_str, ignore_edges=3):
    ref_pos = 0 # Relative to slice start
    query_pos = 0
    variants = []
    
    # Ref Seq Len for Edge Check
    ref_len = len(ref_seq_str)
    
    for op, length in cigar_tuples:
        # 7 = EQ (Match)
        if op == 7: 
            ref_pos += length
            query_pos += length
            
        # 8 = X (Mismatch / SNP)
        elif op == 8:
            for k in range(length):
                # Edge Check
                if ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges):
                    pass # Skip
                else:
                    r_base = ref_seq_str[ref_pos]
                    q_base = query_seq[query_pos + k]
                    variants.append(f"{r_base}{ref_pos + 1}{q_base}")
                
                ref_pos += 1
            query_pos += length
            
        # 1 = I (Insertion)
        elif op == 1:
            if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                inserted_bases = query_seq[query_pos : query_pos + length]
                variants.append(f"ins{ref_pos}{inserted_bases}")
            query_pos += length
            
        # 2 = D (Deletion)
        elif op == 2:
            for k in range(length):
                if not (ref_pos < ignore_edges or ref_pos >= (ref_len - ignore_edges)):
                    variants.append(f"del{ref_pos + 1}")
                ref_pos += 1

        # 4 = S (Soft Clip)
        elif op == 4:
            query_pos += length
            
        # 0 = M (Mixed - should not occur with --eqx)
        elif op == 0:
            ref_pos += length
            query_pos += length

    if not variants:
        return "Ref"
        
    return ",".join(variants)


# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(description="QONTAS: Hybrid VSEARCH+Minimap Amplicon Pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    parser.add_argument("-b", "--base", required=True, help="Output Basename")
    parser.add_argument("-o", "--outdir", default="Qontas_Hybrid_Out", help="Output Directory")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--region", required=True, help="BED file")
    parser.add_argument("--mincount", type=int, default=10, help="Minimum count for a unique sequence")
    parser.add_argument("--ignore_edges", type=int, default=3, help="Ignore variants at edges")
    
    # Global length defaults
    parser.add_argument("--default_minlen", type=int, default=50)
    parser.add_argument("--default_maxlen", type=int, default=10000)

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    regions = parse_bed_regions(args.region, args.default_minlen, args.default_maxlen)
    print(f"[*] Processing {len(regions)} regions.")

    # 1. Map Everything (Global)
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    if not os.path.exists(full_bam):
        print("--- Mapping reads (Global) ---")
        map_cmd = (f"seqkit seq -m 20 {args.input} | " 
                   f"minimap2 -a -x map-ont -t {args.threads} {args.ref} - | "
                   f"samtools view -u -F 2308 | "
                   f"samtools sort -@ {args.threads} -o {full_bam}")
        run_cmd(map_cmd, "Mapping")
        run_cmd(f"samtools index {full_bam}", "Indexing")
    else:
        print("    Using existing BAM.")

    # 2. Loop Regions
    import pysam # Import here to ensure availability
    bam = pysam.AlignmentFile(full_bam, "rb")

    for reg in regions:
        chrom = reg['chrom']
        start = reg['start']
        end = reg['end']
        name = reg['name']
        min_l = reg['min_len']
        max_l = reg['max_len']
        
        print(f"--- Processing {name} ---")
        prefix = f"{args.base}_{name}"
        region_dir = os.path.join(args.outdir, prefix)
        
        # A. Extract Reads to FASTA
        reg_fasta = os.path.join(args.outdir, f"{prefix}.fasta")
        
        try:
            iter_reads = bam.fetch(contig=chrom, start=start, stop=end)
        except ValueError:
            print("    [!] Chrom not in BAM.")
            continue

        processed_count = 0
        with open(reg_fasta, 'w') as f_out:
            for read in iter_reads:
                if read.is_unmapped: continue
                # Coverage & Length Check
                if read.reference_start > start or read.reference_end < end: continue
                if read.query_length < min_l or read.query_length > max_l: continue
                
                # Write to FASTA
                # Note: No reverse complement needed, we just want the raw sequence 
                # to feed into vsearch. Vsearch handles orientation or treats as is.
                # Actually, standard practice is to orient them.
                # Pysam .query_sequence IS oriented to match the reference (forward).
                f_out.write(f">{read.query_name}\n{read.query_sequence}\n")
                processed_count += 1
        
        if processed_count == 0:
            print("    No reads found.")
            continue

        # B. Dereplicate (VSEARCH) - The Denoising Step!
        reg_derep = os.path.join(args.outdir, f"{prefix}_derep.fasta")
        reg_uc = os.path.join(args.outdir, f"{prefix}_derep.txt")
        
        # minuniquesize filters out singleton errors
        run_cmd(f"vsearch --derep_fulllength {reg_fasta} --output {reg_derep} --uc {reg_uc} --minuniquesize {args.mincount} --sizeout",
                f"Dereplicating {name}")

        # C. Map Unique Sequences to Ref Slice (Minimap2 --eqx)
        # Create Ref Slice
        ref_slice_fa = os.path.join(args.outdir, f"{prefix}_ref_slice.fasta")
        if not create_ref_slice_fasta(args.ref, chrom, start, end, ref_slice_fa): continue
        
        unique_bam = os.path.join(args.outdir, f"{prefix}_unique.bam")
        map_u_cmd = (f"minimap2 -a --eqx -x map-ont -t {args.threads} {ref_slice_fa} {reg_derep} | "
                     f"samtools view -b -o {unique_bam}")
        run_cmd(map_u_cmd, "Mapping Unique Sequences")
        
        # D. Parse CIGARs & Build Table
        # Load Ref Slice String
        ref_slice_str = str(next(SeqIO.parse(ref_slice_fa, "fasta")).seq)
        
        # Parse Unique BAM
        ubam = pysam.AlignmentFile(unique_bam, "rb")
        variant_map = {} # ClusterID -> VariantString
        
        for read in ubam.fetch():
            if read.is_unmapped: 
                variant_map[read.query_name] = "Unaligned"
                continue
            
            # Note: The read name in derep fasta might have ";size=N"
            # Vsearch usually outputs ">ID;size=N" in fasta, but ">ID" in BAM query_name?
            # Pysam usually splits at space. Vsearch output usually has no space.
            # We match strictly.
            
            # Clean ID: vsearch might append info, but usually Minimap takes the whole string.
            # If vsearch output >uuid;size=10, minimap BAM query_name is "uuid;size=10"
            clean_id = read.query_name.split(';')[0] # Remove size tag if needed for mapping back to UC?
            # Actually, let's keep the map key consistent with what's in the UC file.
            
            var_str = parse_cigar_for_variants(read.cigartuples, read.query_sequence, ref_slice_str, args.ignore_edges)
            variant_map[read.query_name] = var_str
            # Also map the cleaned ID just in case
            variant_map[clean_id] = var_str
            
        ubam.close()
        
        # E. Process UC File
        try:
            df = pd.read_csv(reg_uc, sep="\t", header=None, dtype=str, on_bad_lines='skip')
        except: continue
        
        if df.empty: continue
        
        # Filter Clusters
        clusters = df[df[0] == "C"][[8, 2]].copy()
        clusters.columns = ["Cluster", "Count"]
        clusters["Count"] = pd.to_numeric(clusters["Count"])
        
        # Map Variants
        # UC Cluster col usually matches the fasta ID.
        clusters["Variant"] = clusters["Cluster"].map(variant_map).fillna("Unknown")
        
        # Aggregate by Variant (Merge identical variants that were split by VSEARCH due to length diffs etc)
        final_df = clusters.groupby("Variant")["Count"].sum().reset_index()
        final_df["RelAbund"] = (final_df["Count"] / final_df["Count"].sum()) * 100
        final_df = final_df.sort_values(by="Count", ascending=False)
        
        # Save
        out_tsv = os.path.join(args.outdir, f"{prefix}_variants.tsv")
        final_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"    Saved: {out_tsv}")
        
        # Cleanup temps
        for f in [ref_slice_fa, unique_bam, reg_derep, reg_uc, reg_fasta]:
            if os.path.exists(f): os.remove(f)

    bam.close()
    print("\nSUCCESS.")

if __name__ == "__main__":
    main()
    