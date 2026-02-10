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
from Bio.Align import PairwiseAligner

# --- Helper for Shell Commands ---
def run_cmd(cmd, description):
    print(f"--- {description} ---")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[!] ERROR: Command failed during: {description}")
        print(f"FAILED COMMAND: {cmd}")
        sys.exit(1)

# --- Helper: Parse BED File (Multi-Region Support) ---
def parse_bed_regions(bed_file):
    regions = []
    try:
        with open(bed_file, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"): continue
                
                parts = line.split('\t')
                if len(parts) < 3: continue
                
                # Column 1: Chromosome/Contig Name (Must match FASTA header)
                chrom = parts[0].strip()
                start = int(parts[1])
                end = int(parts[2])
                
                # Column 4: Region Name (for output filenames)
                if len(parts) >= 4:
                    name = parts[3].strip()
                else:
                    name = f"region_{i+1:02d}"
                    
                regions.append({'chrom': chrom, 'name': name, 'start': start, 'end': end})
    except Exception as e:
        print(f"Error parsing BED file: {e}")
        sys.exit(1)
    return regions

# --- Helper: Get Reference Slice (Robust Fallback) ---
def get_reference_slice(ref_fasta, chrom_id, start, end):
    # 1. Load all records
    records = list(SeqIO.parse(ref_fasta, "fasta"))
    
    target_record = None
    target_id = chrom_id.strip()

    # 2. Try Exact Match
    for r in records:
        if r.id.strip() == target_id:
            target_record = r
            break
            
    # 3. Fallback: If Ref has only 1 sequence, use it regardless of name mismatch
    if target_record is None and len(records) == 1:
        # print(f"    [Info] BED ID '{target_id}' not found, but Ref has 1 sequence. Using '{records[0].id}' as fallback.")
        target_record = records[0]

    if target_record:
        # Safety check for coordinates
        if end > len(target_record.seq):
            print(f"    [!] Warning: Region end ({end}) > Ref length ({len(target_record.seq)}). Truncating.")
            end = len(target_record.seq)
            
        return str(target_record.seq[start:end])
    else:
        print(f"    [!] ERROR: Chromosome '{target_id}' not found in Ref FASTA. SNPs cannot be called.")
        return None

# --- Helper: Call SNPs & Indels ---
def call_variants(query_seq, ref_seq):
    if query_seq == ref_seq:
        return "Ref"
    
    # Configure Aligner (Favors SNPs over Gaps)
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10  # Heavy penalty to prevent spurious gaps
    aligner.extend_gap_score = -1
    
    alignments = aligner.align(ref_seq, query_seq)
    if not alignments:
        return "Unaligned"
    alignment = alignments[0]
    
    # alignment[0] is Ref (target), alignment[1] is Query
    ref_aligned = alignment[0][0]
    query_aligned = alignment[0][1]
    
    variants = []
    ref_pos_counter = 0 # Tracks position in the original reference slice (1-based)
    
    for r_char, q_char in zip(ref_aligned, query_aligned):
        
        # Case 1: Insertion (Gap in Ref)
        if r_char == '-':
            # Insertions happen *after* the current ref base
            variants.append(f"ins{ref_pos_counter}{q_char}")
            continue
            
        # We only advance ref counter if we have a real ref base
        ref_pos_counter += 1
        
        # Case 2: Deletion (Gap in Query)
        if q_char == '-':
            variants.append(f"del{ref_pos_counter}")
            continue
            
        # Case 3: Substitution (SNP)
        if r_char != q_char:
            variants.append(f"{r_char}{ref_pos_counter}{q_char}")

    if not variants:
        return "Ref" 
        
    return ",".join(variants)

# --- 1. Filter FASTQ (Using SeqKit) ---
def filter_fastq_seqkit(input_file, output_file, min_len, max_len, threads):
    print(f"--- Filtering FASTQ (Min: {min_len}, Max: {max_len}) ---")
    cmd = (f"seqkit seq -j {threads} -m {min_len} -M {max_len} {input_file} -o {output_file}")
    run_cmd(cmd, "SeqKit Filtering")

# --- 2. Extract & Trim Reads (BAM Index Optimized) ---
def extract_reads_for_region(input_bam, output_fasta, region_name, chrom, start, end):
    print(f"--- Extracting Region: {region_name} ({chrom}:{start}-{end}) ---")
    
    bam = pysam.AlignmentFile(input_bam, "rb")
    seq_records = []
    processed = 0
    skipped = 0

    # Fetch ONLY reads overlapping this region (Requires Indexed BAM)
    try:
        iter_reads = bam.fetch(contig=chrom, start=start, stop=end)
    except ValueError:
        # Fallback if chrom name in BAM differs from BED
        iter_reads = bam.fetch()

    for read in iter_reads:
        if read.is_unmapped or read.query_sequence is None:
            continue

        # Skip if read doesn't FULLY cover region (Strict Mode)
        if read.reference_start > start or read.reference_end < end:
            skipped += 1
            continue

        q_start = None
        q_end = None
        
        # Robust CIGAR mapping
        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        for q_pos, r_pos in aligned_pairs:
            if q_start is None and r_pos >= start:
                q_start = q_pos
            if r_pos <= end:
                q_end = q_pos

        if q_start is not None and q_end is not None:
            s_idx = min(q_start, q_end)
            e_idx = max(q_start, q_end)
            
            # Slice
            seq_str = read.query_sequence[s_idx:e_idx]
            seq_obj = Seq(seq_str)
            
            # Reorient
            if read.is_reverse:
                seq_obj = seq_obj.reverse_complement()

            record = SeqRecord(seq_obj, id=read.query_name, description=f"region={region_name}")
            seq_records.append(record)
            processed += 1
        else:
            skipped += 1

    bam.close()
    SeqIO.write(seq_records, output_fasta, "fasta")
    print(f"    Extracted {processed} reads.")

# --- 3. Make Feature Table & Call Variants ---
def make_feature_table(input_uc, derep_fasta, output_prefix, min_count, ref_slice=None):
    print(f"--- Generating Feature Table for {output_prefix} ---")
    
    # Load UC
    try:
        df = pd.read_csv(input_uc, sep="\t", header=None, dtype=str, on_bad_lines='skip')
    except Exception:
        df = pd.read_csv(input_uc, sep="\t", header=None, dtype=str, error_bad_lines=False)

    # Filter Cluster Rows
    if df.empty:
        print("    Warning: UC file is empty. No features found.")
        return

    clusters = df[df[0] == "C"][[8, 2]].copy()
    clusters.columns = ["Cluster", "Count"]

    # Convert Count to Numeric
    clusters["Count"] = pd.to_numeric(clusters["Count"])
    
    # Filter by Count
    clusters = clusters[clusters["Count"] >= min_count]
    
    # Rel Abundance
    total = clusters["Count"].sum()
    if total > 0:
        clusters["RelAbund"] = 100 * clusters["Count"] / total
    else:
        clusters["RelAbund"] = 0.0

    # Call Variants
    if ref_slice:
        print("    Calling variants (SNPs/Indels)...")
        # Load sequences into dict for fast lookup
        seq_dict = {r.id: str(r.seq) for r in SeqIO.parse(derep_fasta, "fasta")}
        variants = []
        for cid in clusters["Cluster"]:
            s = seq_dict.get(str(cid))
            if s:
                variants.append(call_variants(s, ref_slice))
            else:
                variants.append("Unknown")
        clusters["Variant"] = variants
    else:
        clusters["Variant"] = "NA"
    
    # Write TSV (Cluster, Variant, Count, RelAbund)
    tsv_out = f"{output_prefix}_table.tsv"
    
    # Reorder columns
    cols = ["Cluster", "Variant", "Count", "RelAbund"] if "Variant" in clusters.columns else ["Cluster", "Count", "RelAbund"]
    clusters[cols].to_csv(tsv_out, sep="\t", index=False)
    print(f"    Table written to {tsv_out}")

    # Write Unique Sequences
    valid_ids = set(clusters['Cluster'].astype(str))
    seqs_out = f"{output_prefix}_unique_sequences.fasta"
    count = 0
    with open(seqs_out, "w") as out_f:
        for r in SeqIO.parse(derep_fasta, "fasta"):
            if r.id in valid_ids:
                if ref_slice and "Variant" in clusters.columns:
                    try:
                        var = clusters.loc[clusters['Cluster'] == r.id, 'Variant'].values[0]
                        r.description += f" variant={var}"
                    except IndexError:
                        pass
                SeqIO.write(r, out_f, "fasta")
                count += 1
    print(f"    Unique sequences: {seqs_out} ({count})")


# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(description="QONTAS: Quantification of ONT Amplicon Sequences")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file")
    parser.add_argument("-b", "--base", required=True, help="Output Basename")
    parser.add_argument("-o", "--outdir", default="Qontas_Out", help="Output Directory (default: Qontas_Out)")
    parser.add_argument("--minlen", type=int, default=600, help="FASTQ reads shorter than this will be discarded (default: 600)")
    parser.add_argument("--maxlen", type=int, default=650, help="FASTQ reads longer than this will be discarded (default: 650)")
    parser.add_argument("--mincount", type=int, default=10, help="Unique sequences must be present this number of times to be reported (default: 10)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--region", required=True, help="BED file (Chrom Start End [Name])")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Parse Regions
    regions = parse_bed_regions(args.region)
    print(f"[*] Found {len(regions)} regions in BED file.")

    # 2. Filter FASTQ (Once)
    filt_fastq = os.path.join(args.outdir, f"{args.base}_filtered.fastq.gz")
    if not os.path.exists(filt_fastq):
        filter_fastq_seqkit(args.input, filt_fastq, args.minlen, args.maxlen, args.threads)
    else:
        print("    Filtered FASTQ exists, skipping step.")

    # 3. Map (Once)
    full_bam = os.path.join(args.outdir, f"{args.base}_full.bam")
    if not os.path.exists(full_bam):
        map_cmd = (f"minimap2 -a -x map-ont -t {args.threads} {args.ref} {filt_fastq} | "
                   f"samtools view -u -F 2308 | "
                   f"samtools sort -@ {args.threads} -o {full_bam}")
        run_cmd(map_cmd, "Mapping")
        run_cmd(f"samtools index {full_bam}", "Indexing")
    else:
        print("    BAM file exists, skipping mapping.")

    # 4. Loop per Region
    for reg in regions:
        r_chrom = reg['chrom']
        r_name = reg['name']
        r_start = reg['start']
        r_end = reg['end']
        
        # ALWAYS suffix region name, even if there is only 1
        prefix = f"{args.base}_{r_name}"
        region_out_base = os.path.join(args.outdir, prefix)
        
        # A. Extract
        reg_fasta = f"{region_out_base}.fasta"
        # Pass CHROM here so we can use bam.fetch() for speed
        extract_reads_for_region(full_bam, reg_fasta, r_name, r_chrom, r_start, r_end)
        
        if os.path.getsize(reg_fasta) == 0:
            print(f"    [!] No reads found for {r_name}. Skipping.")
            continue

        # B. Dereplicate
        reg_derep = f"{region_out_base}_derep.fasta"
        reg_uc = f"{region_out_base}_derep.txt"
        # Removed --threads argument here (VSEARCH derep is single-threaded)
        run_cmd(f"vsearch --derep_fulllength {reg_fasta} --output {reg_derep} --uc {reg_uc}",
                f"Dereplicating {r_name}")

        # C. Feature Table & Variants
        # Pass CHROM here so we find the right ref sequence
        ref_slice = get_reference_slice(args.ref, r_chrom, r_start, r_end)
        make_feature_table(reg_uc, reg_derep, region_out_base, args.mincount, ref_slice)

    # Cleanup
    print("\n--- Cleanup ---")
    if os.path.exists(full_bam): os.remove(full_bam)
    if os.path.exists(full_bam+".bai"): os.remove(full_bam+".bai")
    
    print(f"\nSUCCESS. Results in {args.outdir}")

if __name__ == "__main__":
    main()
    