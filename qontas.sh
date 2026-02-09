#!/usr/bin/env bash

echo ''
INPUTREADS=$1
REFERENCE=$2
BASENAME=$3
OUTDIR=$4
FILT_MINLEN=$5
FILT_MAXLEN=$6
FEAT_MINCOUNT=$7
FEAT_MINABUND=$8
THREADS=$9
REGION=${10}

echo 'Input FASTQ is' $1
echo 'Reference FASTA is' $2
echo 'Output basename is' $3
echo 'Output directory is' $4
echo 'FASTQ reads shorter than' $5 'bp will be discarded'
echo 'FASTQ reads longer than' $6 'bp will be discarded'
echo 'Sequences observed fewer than' $7 'times per sample will be discarded BEFORE calculating relative abundance'
echo 'Sequences with per-sample relative abundance below' $8 'will be discarded'
echo 'Using' $9 'threads for minimap2 mapping'

if [ -n "$REGION" ]; then
    echo "Extracting regions defined in: $REGION"
else
    echo "No BED file provided. Processing all mapped reads."
fi

echo ''

mkdir -p "$OUTDIR"

filter_fastq.py \
    --input_file "$INPUTREADS" \
    --output_file "$OUTDIR"/"$BASENAME"_filtered.fastq.gz \
    --min_length "$FILT_MINLEN" \
    --max_length "$FILT_MAXLEN"

minimap2 \
    -a \
    -x map-ont \
    "$REFERENCE" \
    "$OUTDIR"/"$BASENAME"_filtered.fastq.gz | \
    samtools view -u -F 2308 | \
    samtools sort -@ "$THREADS" -O BAM -o "$OUTDIR"/"$BASENAME"_full.bam

samtools index "$OUTDIR"/"$BASENAME"_full.bam

if [ -z "$REGION" ]; then

    # not region filtered, use the full BAM
    mv "$OUTDIR"/"$BASENAME"_full.bam "$OUTDIR"/"$BASENAME".bam
    mv "$OUTDIR"/"$BASENAME"_full.bam.bai "$OUTDIR"/"$BASENAME".bam.bai

else

    # extract specified region from BAM
    samtools view -b -L "$REGION" "$OUTDIR"/"$BASENAME"_full.bam > "$OUTDIR"/"$BASENAME".bam
    samtools index "$OUTDIR"/"$BASENAME".bam

fi

samtools index "$OUTDIR"/"$BASENAME".bam

write_forward_reads.py \
    --input_bam "$OUTDIR"/"$BASENAME".bam \
    --output_fasta "$OUTDIR"/"$BASENAME".fasta

vsearch \
    --derep_fulllength "$OUTDIR"/"$BASENAME".fasta \
    --output "$OUTDIR"/"$BASENAME"_derep.fasta \
    --uc "$OUTDIR"/"$BASENAME"_derep.txt

make_feature_table.py \
    -i "$OUTDIR"/"$BASENAME"_derep.txt \
    -o "$OUTDIR"/"$BASENAME" \
    --min_count "$FEAT_MINCOUNT" \
    --min_abundance "$FEAT_MINABUND"

seqkit grep -f "$OUTDIR"/"$BASENAME"_readnames.txt \
    "$OUTDIR"/"$BASENAME".fasta > "$OUTDIR"/"$BASENAME"_sequences.fasta

rm -f "$OUTDIR"/"$BASENAME"_filtered.fastq.gz \
    "$OUTDIR"/"$BASENAME"_full.bam "$OUTDIR"/"$BASENAME"_full.bam.bai \
    "$OUTDIR"/"$BASENAME".bam "$OUTDIR"/"$BASENAME".bam.bai \
    "$OUTDIR"/"$BASENAME".fasta \
    "$OUTDIR"/"$BASENAME"_derep.fasta "$OUTDIR"/"$BASENAME"_derep.txt \
    "$OUTDIR"/"$BASENAME"_readnames.txt
