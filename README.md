# Qontas  
### :airplane: Quantification of Oxford Nanopore Techonologies Amplicon Sequences :airplane:   

[![DOI](https://zenodo.org/badge/652407080.svg)](https://zenodo.org/badge/latestdoi/652407080)  

:construction: under construction :construction:  

The aim of this tool is to take amplicon sequencing reads and generate a table of all variants identified, along with the counts and relative abundance abundances.  

Input:
1. Sequencing reads in FASTQ format (can be gzipped)
2. Reference sequence of the target region in FASTA format.  

Output:   
1. Three column TSV file listing, for unique sequence, the name, count, and relative abundance (`basename_clusters.txt`)
2. MultiFASTA file containing the sequence of each unique sequence named in the TSV output file (`basename_sequences.fasta`)

## SETUP  

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

## USAGE  

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

## TO DO
- [ ] automate pipeline in (better) single executable
- [ ] give option to retain or detele temp files
- [ ] give option to specify output directory
- [ ] modify to accept a list of input FASTQ files (TSV format) and output a single merged feature table  
   Will need to modify to generate md5 read names so that they are groupable between samples  
- [ ] include flag to modify minimap -x flag allowing PacBio (`-x map-pb`) or Illumina (`-x sr`) reads