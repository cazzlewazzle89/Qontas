#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import pysam
import pandas as pd
from collections import Counter
from Bio import SeqIO

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
                
                # Custom Lengths per Region (Cols 5 & 6)
                # Fallback to defaults if not specified
                if len(parts) >= 6:
                    r_min = int(parts[4])
                    r_max = int(parts[5])
                else:
                    r_min = default_min
                    r_max = default_max

                regions.append({
                    'chrom': chrom, 
                    'name': name, 
                    'start': start, 
                    'end': end,
                    'min_len': r_min,
                    'max_len': r_max
                })
    except Exception as e:
        print(f"Error parsing BED file: {e}")
        sys.exit(1)
    return regions

# --- Parse CIGAR for a specific region ---
def get_variants_in_region(read, region_start, region_end, ref_seq_str, ignore_edges=3):
    ref_pos = read.reference_start
    query_pos = 0
    variants = []
    
    cigar_tuples = read.cigartuples
    if not cigar_tuples: return "Unaligned"

    if ref_pos > region_start: return None 
    
    for op, length in cigar_tuples:
        if ref_pos >= region_end: break

        # 7 = EQ (Match)
        if op == 7: 
            ref_pos += length
            query_pos += length
            
        # 8 = X (Mismatch / SNP)
        elif op == 8:
            for k in range(length):
                current_r_pos = ref_pos + k
                if region_start <= current_r_pos < region_end:
                    dist_from_start = current_r_pos - region_start
                    dist_from_end = region_end - current_r_pos
                    if not (dist_from_start < ignore_edges or dist_from_end <= ignore_edges):
                        r_base = ref_seq_str[current_r_pos]
                        q_base = read.query_sequence[query_pos + k]
                        rel_pos = current_r_pos - region_start + 1
                        variants.append(f"{r_base}{rel_pos}{q_base}")
            ref_pos += length
            query_pos += length

        # 1 = I (Insertion)
        elif op == 1:
            if region_start <= ref_pos < region_end:
                dist_from_start = ref_pos - region_start
                dist_from_end = region_end - ref_pos
                if not (dist_from_start < ignore_edges or dist_from_end <= ignore_edges):
                    inserted_bases = read.query_sequence[query_pos : query_pos + length]
                    rel_pos = ref_pos - region_start
                    variants.append(f"ins{rel_pos}{inserted_bases}")
            query_pos += length
            
        # 2 = D (Deletion)
        elif op == 2:
            for k in range(length):
                current_r_pos = ref_pos + k
                if region_start <= current_r_pos < region_end:
                    dist_from_start = current_r_pos - region_start
                    dist_from_end = region_end - current_r_pos
                    if not (dist_from_start < ignore_edges or dist_from_end <= ignore_edges):
                        rel_pos = current_r_pos - region_start + 1
                        variants.append(f"del{rel_pos}")
            ref_pos += length

        # 4 = S (Soft Clip)
        elif op == 4:
            query_pos += length
            
        # 0 = M (Mixed)
        elif op == 0:
            ref_pos += length
            query_pos += length

    if not variants: return "Ref"
    return ",".join(variants)


# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(description="QONTAS: Multiplexed Amplicon Quantification")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    parser.add_argument("-b", "--base", required=True, help="Output Basename")
    parser.add_argument("-o", "--outdir", default="Qontas_Direct_Out", help="Output Directory")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--region", required=True, help="BED file (Chrom Start End Name [MinLen MaxLen])")
    parser.add_argument("--ignore_edges", type=int, default=3, help="Ignore variants at edges")
    
    # Global defaults (used if BED columns 5/6 are missing)
    parser.add_argument("--default_minlen", type=int, default=50, help="Default min read length if not in BED")
    parser.add_argument("--default_maxlen", type=int, default=10000, help="Default max read length if not in BED")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Parse Regions & Lengths
    regions = parse_bed_regions(args.region, args.default_minlen, args.default_maxlen)
    print(f"[*] Processing {len(regions)} regions.")
    for r in regions:
        print(f"    - {r['name']}: {r['chrom']}:{r['start']}-{r['end']} (Length Filter: {r['min_len']}-{r['max_len']} bp)")

    ref_genome = {r.id: str(r.seq) for r in SeqIO.parse(args.ref, "fasta")}

    # 2. Map Everything (No strict pre-filtering)
    # We map ALL reads. We only filter extremely short junk (<20bp) to save time/space.
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    
    if not os.path.exists(full_bam):
        print("--- Mapping reads (Global) ---")
        # Note: Added -m 20 to seqkit just to clear empty/garbage reads
        map_cmd = (f"seqkit seq -m 20 {args.input} | " 
                   f"minimap2 -a --eqx -x map-ont -t {args.threads} {args.ref} - | "
                   f"samtools view -u -F 2308 | "
                   f"samtools sort -@ {args.threads} -o {full_bam}")
        run_cmd(map_cmd, "Mapping")
        run_cmd(f"samtools index {full_bam}", "Indexing")
    else:
        print("    Using existing BAM.")

    # 3. Quantify Regions
    bam = pysam.AlignmentFile(full_bam, "rb")

    for reg in regions:
        chrom = reg['chrom']
        start = reg['start']
        end = reg['end']
        name = reg['name']
        
        # Region-specific length thresholds
        min_l = reg['min_len']
        max_l = reg['max_len']
        
        print(f"--- Quantifying {name} ---")
        
        if chrom not in ref_genome:
            print(f"    [!] Error: Chrom '{chrom}' not in reference FASTA. Skipping.")
            continue
        ref_seq_str = ref_genome[chrom]

        variant_counts = Counter()
        read_processed = 0
        read_skipped_cov = 0
        read_skipped_len = 0
        
        try:
            iter_reads = bam.fetch(contig=chrom, start=start, stop=end)
        except ValueError:
            print("    [!] Chromosome not found in BAM.")
            continue

        for read in iter_reads:
            if read.is_unmapped: continue
            
            # A. Check Coverage
            if read.reference_start > start or read.reference_end < end:
                read_skipped_cov += 1
                continue
            
            # B. Check Length (Region Specific!)
            # query_length returns the length of the sequence in the BAM
            r_len = read.query_length
            if r_len < min_l or r_len > max_l:
                read_skipped_len += 1
                continue
            
            # C. Call Variants
            var_str = get_variants_in_region(read, start, end, ref_seq_str, args.ignore_edges)
            
            if var_str:
                variant_counts[var_str] += 1
                read_processed += 1
        
        print(f"    Accepted: {read_processed}")
        print(f"    Filtered (Coverage): {read_skipped_cov}")
        print(f"    Filtered (Length {min_l}-{max_l}): {read_skipped_len}")

        # 4. Output Table
        if not variant_counts:
            print("    No reads passed filters.")
            continue
            
        df = pd.DataFrame.from_dict(variant_counts, orient='index', columns=['Count'])
        df.index.name = 'Variant'
        df.reset_index(inplace=True)
        
        total = df['Count'].sum()
        df['RelAbund'] = (df['Count'] / total) * 100
        df = df.sort_values(by='Count', ascending=False)
        
        out_file = os.path.join(args.outdir, f"{args.base}_{name}_variants.tsv")
        df.to_csv(out_file, sep="\t", index=False)
        print(f"    Saved: {out_file}")

    bam.close()
    print("\nSUCCESS.")

if __name__ == "__main__":
    main()
