#!/usr/bin/env bash

INPUTREADS='/home/aminichiello/testingAmpliconSequences/DMG2303177.fastq.gz'
BASENAME='DMG2303177'
REFERENCE='/home/aminichiello/testingAmpliconSequences/rpoBreference.fasta'
FILT_MINLEN=600
FILT_MAXLEN=650
FEAT_MINCOUNT=2
FEAT_MINABUND=0.1
THREADS=10

filter_fastq.py \
    --input_file "$INPUTREADS" \
    --output_file "$BASENAME"_filtered.fastq.gz \
    --min_length "$FILT_MINLEN" \
    --max_length "$FILT_MAXLEN"

minimap2 \
    -a \
    -x map-ont \
    "$REFERENCE" \
    "$BASENAME"_filtered.fastq.gz | samtools \
        sort | samtools \
            view \
            -b \
            -F 4 > "$BASENAME".bam

samtools index "$BASENAME".bam

write_forward_reads.py \
    --bam_file "$BASENAME".bam \
    --output_fasta "$BASENAME".fasta

vsearch \
    --derep_fulllength "$BASENAME".fasta \
    --output "$BASENAME"_derep.fasta \
    --threads "$THREADS" \
    --uc "$BASENAME"_derep.txt

make_feature_table.py \
    -i "$BASENAME"_derep.txt \
    -o "$BASENAME"_clusters.txt \
    --min_count "$FEAT_MINCOUNT" \
    --min_abundance "$FEAT_MINABUND"

cut -f 1 "$BASENAME"_clusters.txt | sed '1d' > "$BASENAME"_readnames.txt

seqkit grep -f "$BASENAME"_readnames.txt "$BASENAME".fasta > "$BASENAME"_sequences.fasta