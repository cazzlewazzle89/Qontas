# Qontas  
### Quantification of Oxford Nanopore Techonologies Amplicon Sequences   

Disclaimer: :construction: under construction :construction:  
 
Create conda environment
```bash
conda create -n ont_amplicon -c bioconda -c conda-forge minimap2 samtools=1.11 seqiolib seqkit vsearch pandas -y
```

Clone repo and make scripts executable
 ```bash
git clone https://github.com/cazzlewazzle89/Qontas.git

chmod +x Qontas/*
```

Add directory (eg. `/home/cwwalsh/Software/Qontas`) to your path  
Handy guide [here](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/) 

Run pipeline using the following commands - currently working on automating this
```bash
# take raw reads and remove those too long/short (likely off-target amplicons)
# set appropriate --min_length and --max_length arguments based on your expected amplicon size
# doesn't need to be too strict, just filter out obvious artifacts
filter_fastq.py --input_file sample.fastq.gz --output_file sample_filtered.fastq.gz --min_length 600 --max_length 650

# align reads to a reference sequence to remove any further off-target amplicons
minimap2 -a -x map-ont reference.fa sample_filtered.fastq.gz | samtools sort | samtools view -b -F 4 > sample.bam

# index BAM file
samtools index sample.bam

# orient reads in "correct" forward orientation
write_forward_reads.py --bam_file sample.bam --output_fasta sample.fasta

# identify unique reads and count their frequency
vsearch --derep_fulllength sample.fasta --output sample_derep.fasta --threads 10 --uc sample_derep.txt

# convert this information to a long-form (tidyverse-friendly) feature table
# can modify min_count and min_abundance flags
# see --help for more information 
# will output a 3 column TSV, listing: sequenceID, count, and relative abundance for each unique seqeunce
make_feature_table.py -i sample_derep.txt -o sample_clusters.txt --min_count 2 --min_abundance 0.1

# print the squenceIDs to a new file
cut -f 1 sample_clusters.txt | sed '1d' > sample_readnames.txt

# extact these sequences in FASTA format
seqkit grep -f sample_readnames.txt sample.fasta > sample_sequences.fasta
```

## TO DO
- [ ] automate pipeline in single executable
- [ ] modify to accept a list of input FASTA files (TSV format) and output a single merged feature table
- [ ] include flag modifying minimap -x flag allow PacBio (`-x map-pb`) or Illumina (`-x sr`) reads