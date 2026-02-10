#!/usr/bin/env python3
import argparse, subprocess, os, sys, pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd, desc):
    print(f"--- {desc} ---")
    subprocess.run(cmd, shell=True, check=True)

def parse_bed(bed_file):
    regs = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.strip().split('\t')
            regs.append({'chrom': p[0], 'start': int(p[1]), 'end': int(p[2]), 'name': p[3]})
    return regs

def get_variants(read, ref_seq, ignore=5):
    """Parses =/X CIGAR for SNVs, ignoring N bases at the ends."""
    if not read.cigartuples: return "Unaligned"
    ref_pos, q_pos, vars = 0, 0, []
    ref_len = len(ref_seq)
    for op, length in read.cigartuples:
        if op == 7: # Match (=)
            ref_pos += length; q_pos += length
        elif op == 8: # Mismatch (X)
            for k in range(length):
                if ignore <= ref_pos < (ref_len - ignore):
                    vars.append(f"{ref_seq[ref_pos]}{ref_pos+1}{read.query_sequence[q_pos+k]}")
                ref_pos += 1; q_pos += 1
        elif op == 1: # Ins
            if ignore <= ref_pos < (ref_len - ignore):
                vars.append(f"ins{ref_pos}{read.query_sequence[q_pos:q_pos+length]}")
            q_pos += length
        elif op == 2: # Del
            for k in range(length):
                if ignore <= ref_pos < (ref_len - ignore):
                    vars.append(f"del{ref_pos+1}")
                ref_pos += 1
        else: # S, H, M
            ref_pos += length if op in [0, 2, 3, 7, 8] else 0
            q_pos += length if op in [0, 1, 4, 7, 8] else 0
    return ",".join(vars) if vars else "Ref"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-b", "--bed", required=True)
    parser.add_argument("-o", "--outdir", default="Qontas_v2")
    parser.add_argument("--mincount", type=int, default=10)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Align all reads once
    bam = os.path.join(args.outdir, "initial.bam")
    run_cmd(f"minimap2 -ax map-ont {args.ref} {args.input} | samtools sort -o {bam}", "Global Alignment")
    run_cmd(f"samtools index {bam}", "Indexing")

    regions = parse_bed(args.bed)
    ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))

    for reg in regions:
        print(f"\nProcessing {reg['name']}...")
        clipped_fa = os.path.join(args.outdir, f"{reg['name']}_clipped.fa")
        
        # 2. Extract and Hard Clip
        with pysam.AlignmentFile(bam, "rb") as sam, open(clipped_fa, "w") as f:
            for r in sam.fetch(reg['chrom'], reg['start'], reg['end']):
                if r.reference_start <= reg['start'] and r.reference_end >= reg['end']:
                    # Use get_aligned_pairs to find the specific query bases for the BED start/end
                    pairs = r.get_aligned_pairs(matches_only=True)
                    q_s = next((q for q, r_p in pairs if r_p >= reg['start']), None)
                    q_e = next((q for q, r_p in reversed(pairs) if r_p <= reg['end']), None)
                    if q_s is not None and q_e is not None:
                        s, e = (q_s, q_e) if q_s < q_e else (q_e, q_s)
                        seq = Seq(r.query_sequence[s:e+1])
                        if r.is_reverse: seq = seq.reverse_complement()
                        f.write(f">{r.query_name}\n{str(seq)}\n")

        # 3. Denoise with VSEARCH
        derep_fa = os.path.join(args.outdir, f"{reg['name']}_derep.fa")
        uc_file = os.path.join(args.outdir, f"{reg['name']}_derep.txt")
        run_cmd(f"vsearch --derep_fulllength {clipped_fa} --output {derep_fa} --uc {uc_file} --minuniquesize {args.mincount} --sizeout", "Denoising")

        # 4. Map Unique Seqs back to Ref Slice
        ref_slice = ref_dict[reg['chrom']].seq[reg['start']:reg['end']]
        slice_fa = os.path.join(args.outdir, f"{reg['name']}_slice.fa")
        SeqIO.write(SeqRecord(ref_slice, id="slice"), slice_fa, "fasta")
        
        ubam = os.path.join(args.outdir, f"{reg['name']}_unique.bam")
        run_cmd(f"minimap2 -ax map-ont --eqx {slice_fa} {derep_fa} | samtools sort -o {ubam}", "Final Naming Alignment")
        
        # 5. Final Table
        var_map = {}
        with pysam.AlignmentFile(ubam, "rb") as u:
            for r in u:
                if not (r.is_secondary or r.is_supplementary):
                    var_map[r.query_name.split(';')[0]] = get_variants(r, str(ref_slice))

        df = pd.read_csv(uc_file, sep="\t", header=None, dtype=str)
        clusters = df[df[0] == 'C'][[8, 2]].copy()
        clusters.columns = ['ID', 'Count']
        clusters['ID'] = clusters['ID'].str.split(';').str[0]
        clusters['Count'] = pd.to_numeric(clusters['Count'])
        clusters['Variant'] = clusters['ID'].map(var_map).fillna("Unaligned")
        
        final = clusters.groupby("Variant")['Count'].sum().reset_index().sort_values("Count", ascending=False)
        final['RelAbund'] = (final['Count'] / final['Count'].sum()) * 100
        final.to_csv(os.path.join(args.outdir, f"{reg['name']}_results.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()
    