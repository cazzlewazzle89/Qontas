#!/usr/bin/env python3
import argparse, subprocess, os, sys, pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd, desc):
    print(f"--- {desc} ---")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[!] ERROR in {desc}: {e}")
        sys.exit(1)

def parse_bed(bed_file):
    regs = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.strip().split('\t')
            # Columns: 0:Chrom, 1:Start, 2:End, 3:Name, 4:MinLen, 5:MaxLen
            regs.append({
                'chrom': p[0], 
                'start': int(p[1]), 
                'end': int(p[2]), 
                'name': p[3],
                'min_len': int(p[4]) if len(p) >= 5 else 100,
                'max_len': int(p[5]) if len(p) >= 6 else 1500
            })
    return regs

def main():
    parser = argparse.ArgumentParser(description="QONTAS: Amplicon Binning with Filtered Abundance Calculation")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA")
    parser.add_argument("-b", "--bed", required=True, help="BED file with lengths in cols 5 & 6")
    parser.add_argument("-s", "--sample", required=True, help="Sample name for subdirectory")
    parser.add_argument("-o", "--outdir", default="Qontas_Out", help="Parent output directory")
    parser.add_argument("--mincount", type=int, default=10, help="Min reads per cluster")
    
    args = parser.parse_args()
    
    workdir = os.path.join(args.outdir, args.sample)
    os.makedirs(workdir, exist_ok=True)

    # 1. Global Alignment
    bam = os.path.join(workdir, f"{args.sample}_initial.bam")
    bam_index = bam + ".bai"
    if not os.path.exists(bam):
        run_cmd(f"minimap2 -ax map-ont {args.ref} {args.input} | samtools view -u -F 2308 | samtools sort -o {bam}", f"Global Alignment: {args.sample}")
        run_cmd(f"samtools index {bam}", "Indexing")

    regions = parse_bed(args.bed)

    for reg in regions:
        region_name = reg['name']
        print(f"\n>>> Processing Region: {region_name}")
        
        clipped_fa = os.path.join(workdir, f"{region_name}_temp_clipped.fa")
        
        # 2. Extract and Hard Clip
        with pysam.AlignmentFile(bam, "rb") as sam, open(clipped_fa, "w") as f:
            read_count = 0
            try:
                for r in sam.fetch(reg['chrom'], reg['start'], reg['end']):
                    if r.query_sequence is None or r.is_secondary or r.is_supplementary: continue
                    if r.reference_start <= reg['start'] and r.reference_end >= reg['end']:
                        pairs = r.get_aligned_pairs(matches_only=True)
                        q_s = next((q for q, r_p in pairs if r_p >= reg['start']), None)
                        q_e = next((q for q, r_p in reversed(pairs) if r_p <= reg['end']), None)
                        if q_s is not None and q_e is not None:
                            s, e = (q_s, q_e) if q_s < q_e else (q_e, q_s)
                            seq_obj = Seq(r.query_sequence[s:e+1])
                            if r.is_reverse: seq_obj = seq_obj.reverse_complement()
                            f.write(f">{r.query_name}\n{str(seq_obj)}\n")
                            read_count += 1
            except ValueError: continue
            
        if read_count < args.mincount:
            print(f"    [!] Insufficient reads ({read_count}). Skipping.")
            if os.path.exists(clipped_fa): os.remove(clipped_fa)
            continue

        # 3. Clustering
        final_fa = os.path.join(workdir, f"{region_name}.fa")
        final_tsv = os.path.join(workdir, f"{region_name}.tsv")
        uc_file = os.path.join(workdir, f"{region_name}_derep.txt")
        
        vsearch_cmd = (
            f"vsearch --derep_fulllength {clipped_fa} --strand plus "
            f"--minseqlength {reg['min_len']} --maxseqlength {reg['max_len']} "
            f"--output {final_fa} --uc {uc_file} --minuniquesize {args.mincount} --sizeout"
        )
        run_cmd(vsearch_cmd, f"Clustering {region_name}")

        # 4. Table Assembly (Filtered Abundance)
        try:
            df_uc = pd.read_csv(uc_file, sep="\t", header=None, dtype=str)
            # Filter for Cluster 'C' rows
            clusters = df_uc[df_uc[0] == 'C'][[8, 2]].copy()
            clusters.columns = ['ClusterID', 'Count']
            clusters['Count'] = pd.to_numeric(clusters['Count'])
            
            # --- Key Change: Abundance based only on surviving clusters ---
            surviving_total = clusters['Count'].sum()
            clusters['RelAbund'] = (clusters['Count'] / surviving_total) * 100 if surviving_total > 0 else 0
            
            clusters = clusters.sort_values("Count", ascending=False)
            clusters.to_csv(final_tsv, sep="\t", index=False)
            print(f"    [Success] {region_name}.tsv generated from {len(clusters)} clusters.")
        except Exception as e:
            print(f"    [!] Table Error: {e}")

        # Cleanup region-specific temps
        for tmp in [clipped_fa, uc_file]:
            if os.path.exists(tmp): os.remove(tmp)

    # 5. Final Cleanup
    print(f"\n--- Finalizing Sample {args.sample}: Removing initial BAM files ---")
    for f_to_del in [bam, bam_index]:
        if os.path.exists(f_to_del):
            os.remove(f_to_del)

    print(f"DONE. Results for {args.sample} located in {workdir}")

if __name__ == "__main__":
    main()