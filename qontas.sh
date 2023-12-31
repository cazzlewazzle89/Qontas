#!/usr/bin/env bash

echo ''
INPUTREADS=$1
REFERENCE=$2
BASENAME=$3
FILT_MINLEN=$4
FILT_MAXLEN=$5
FEAT_MINCOUNT=$6
FEAT_MINABUND=$7
THREADS=$8

echo 'Input FASTQ is' $1
echo 'Reference FASTA is' $2
echo 'Output basename is' $3
echo 'FASTQ reads shorter than' $4 'bp will be discarded'
echo 'FASTQ reads longer than' $5 'bp will be discarded'
echo 'Sequences observed fewer than' $6 'times per sample will be discarded BEFORE calculating relative abundance'
echo 'Sequences with per-sample relative abundance below' $7 'will be discarded'
echo 'Using' $8 'threads for minimap2 mapping'
echo ''

sleep 10

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
    --input_bam "$BASENAME".bam \
    --output_fasta "$BASENAME".fasta

vsearch \
    --derep_fulllength "$BASENAME".fasta \
    --output "$BASENAME"_derep.fasta \
    --threads "$THREADS" \
    --uc "$BASENAME"_derep.txt

make_feature_table.py \
    -i "$BASENAME"_derep.txt \
    -o "$BASENAME" \
    --min_count "$FEAT_MINCOUNT" \
    --min_abundance "$FEAT_MINABUND"

seqkit grep -f "$BASENAME"_readnames.txt "$BASENAME".fasta > "$BASENAME"_sequences.fasta

rm -f "$BASENAME"_filtered.fastq.gz "$BASENAME".bam "$BASENAME".bam.bai "$BASENAME".fasta "$BASENAME"_derep.fasta "$BASENAME"_derep.txt "$BASENAME"_readnames.txt