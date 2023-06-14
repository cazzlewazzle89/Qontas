# Qontas  
### :airplane: Quantification of Oxford Nanopore Techonologies Amplicon Sequences :airplane:   

[![DOI](https://zenodo.org/badge/652407080.svg)](https://zenodo.org/badge/latestdoi/652407080)  

:construction: under construction :construction:  

The aim of this tool is to take amplicon sequencing reads and generate a table of all variants identified, along with the counts and relative abundance abundances.  
Input:
1. sequence reads in FASTQ format (can be gzipped)
2. reference sequence of the target region in FASTA format.  

Output:   
1. three column TSV file listing, for unique sequence, the name, count, and relative abundance
2. MultiFASTA file containing the sequence of each unique seuqneced named in the TSV output file.

 
Create conda environment
```bash
conda create -n qontas -c bioconda minimap2 phylopandas pysam seqiolib seqkit vsearch -y
```

Clone repo and make scripts executable
 ```bash
git clone https://github.com/cazzlewazzle89/Qontas.git

chmod +x Qontas/*
```

Add directory (eg. `/home/cwwalsh/Software/Qontas`) to your path  
Handy guide [here](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/) 

Current full pipeline is run using the script `qontas.sh` with 8 positional parameters.  
Namely:  
1. Input FASTQ
2. Input reference FASTA
3. Output basename
4. Minimum read length for FASTQ filtering
5. Maximum read length for FASTQ filtering
6. Amount of times a sequence must be observed per sample to be retained for relative abundance calculation (recommended to set this >1 to remove singletons [highly likely to be PCR artifacts or sequencing errors])
7. Mimimum relative abundance as a percentage (eg. 2 or 0.1) for a sequence to be reported 
8. Number of threads to use for minimap2 mapping  

eg. `qontas.sh sample.fastq.gz ref.fa sample 600 650 2 0.1 10`  

`qontas` will print these values to screen and wait 10 seconds before running to give you a chance to cancel if anything is wrong  

Alternatively, you can manually run the pipeline using the following commands
```bash
# take raw reads and remove those too long/short (likely off-target amplicons)
# set appropriate --min_length and --max_length arguments based on your expected amplicon size
# doesn't need to be too strict, just filter out obvious artifacts
filter_fastq.py --input_file sample.fastq.gz --output_file sample_filtered.fastq.gz --min_length 600 --max_length 650

# align reads to a reference sequence to remove any further off-target amplicons
minimap2 -a -x map-ont -t 10 reference.fa sample_filtered.fastq.gz | samtools sort | samtools view -b -F 4 > sample.bam

# index BAM file
samtools index sample.bam

# orient reads in "correct" forward orientation
write_forward_reads.py --bam_file sample.bam --output_fasta sample.fasta

# identify unique reads and count their frequency
vsearch --derep_fulllength sample.fasta --output sample_derep.fasta --uc sample_derep.txt

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
- [ ] automate pipeline in (better) single executable
- [ ] modify to accept a list of input FASTQ files (TSV format) and output a single merged feature table
- [ ] include flag to modify minimap -x flag allowing PacBio (`-x map-pb`) or Illumina (`-x sr`) reads