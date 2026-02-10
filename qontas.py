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
            regs.append({'chrom': p[0], 'start': int(p[1]), 'end': int(p[2]), 'name': p[3]})
    return regs

def main():
    parser = argparse.ArgumentParser(description="QONTAS: Simplified Amplicon Binning & Clustering")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA")
    parser.add_argument("-b", "--bed", required=True, help="BED file")
    parser.add_argument("-s", "--sample", required=True, help="Sample name")
    parser.add_argument("-o", "--outdir", default="Qontas_Simplified")
    parser.add_argument("--mincount", type=int, default=10, help="Min reads per cluster")
    
    args = parser.parse_args()
    workdir = os.path.join(args.outdir, args.sample)
    os.makedirs(workdir, exist_ok=True)

    # 1. Global Alignment
    bam = os.path.join(workdir, f"{args.sample}_initial.bam")
    if not os.path.exists(bam):
        run_cmd(f"minimap2 -ax map-ont {args.ref} {args.input} | samtools view -u -F 2308 | samtools sort -o {bam}", "Global Alignment")
        run_cmd(f"samtools index {bam}", "Indexing")

    regions = parse_bed(args.bed)

    for reg in regions:
        print(f"\n>>> Processing Region: {reg['name']}...")
        clipped_fa = os.path.join(workdir, f"{reg['name']}_clipped.fa")
        
        # 2. Extract and Hard Clip to BED coordinates
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
            continue

        # 3. Denoise/Cluster with VSEARCH
        derep_fa = os.path.join(workdir, f"{reg['name']}_variants.fa")
        uc_file = os.path.join(workdir, f"{reg['name']}_derep.txt")
        # Note: --sizeout is kept here so the FASTA headers contain counts for easy viewing
        run_cmd(f"vsearch --derep_fulllength {clipped_fa} --strand plus --output {derep_fa} --uc {uc_file} --minuniquesize {args.mincount} --sizeout", "Clustering")

        # 4. Generate Final Table from UC file
        try:
            df_uc = pd.read_csv(uc_file, sep="\t", header=None, dtype=str)
            # Column 8 (index 8) is the Cluster ID, Column 2 (index 2) is the Count
            clusters = df_uc[df_uc[0] == 'C'][[8, 2]].copy()
            clusters.columns = ['ClusterID', 'Count']
            clusters['Count'] = pd.to_numeric(clusters['Count'])
            
            # Calculate Abundance
            total_reads = clusters['Count'].sum()
            clusters['RelAbund'] = (clusters['Count'] / total_reads) * 100 if total_reads > 0 else 0
            clusters = clusters.sort_values("Count", ascending=False)
            
            res_file = os.path.join(workdir, f"{args.sample}_{reg['name']}_counts.tsv")
            clusters.to_csv(res_file, sep="\t", index=False)
            
            print(f"    [Success] Found {len(clusters)} unique clusters.")
            print(f"    [Output] Table: {res_file}")
            print(f"    [Output] Sequences: {derep_fa}")
            
        except Exception as e:
            print(f"    [!] Table Error: {e}")

        # 5. Cleanup temporary clipped file
        if os.path.exists(clipped_fa): os.remove(clipped_fa)
        if os.path.exists(uc_file): os.remove(uc_file)

    print(f"\nDONE. Results for {args.sample} are in {workdir}")

if __name__ == "__main__":
    main()