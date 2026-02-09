# Qontas  
### :airplane: Quantification of ONT Amplicon Sequences :airplane:   

[![DOI](https://zenodo.org/badge/652407080.svg)](https://zenodo.org/badge/latestdoi/652407080)  

:construction: under construction :construction:  

The aim of this tool is to take amplicon sequencing reads and generate a table of all variants identified, their counts, and relative abundances.  

Input:
1. Sequencing reads in FASTQ format (can be gzipped)
2. Reference sequence of the target region in FASTA format.

Output:   
`basename_variants.txt` - three column TSV file listing, for unique combination of variants, the variants, count, and relative abundance

## SETUP  

Clone repo and make scripts executable
```bash
git clone https://github.com/cazzlewazzle89/Qontas.git ~/Software/QONTAS

chmod +x ~/Software/QONTAS/*
```

Create conda environment
```bash
mamba create -n qontas bioconda::minimap2 bioconda::samtools bioconda::vsearch bioconda::seqkit bioconda::pysam conda-forge::pandas conda-forge::biopython
```

```

Add directory (eg. `~/Software/QONTAS`) to your path  
Handy guide [here](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/)   

## USAGE  

Current full pipeline is run using the script `qontas.sh` with 8 positional parameters.  
Namely:  
1. Input FASTQ
2. Input reference FASTA
3. Output basename
4. Output directory
5. Minimum read length for FASTQ filtering
6. Maximum read length for FASTQ filtering
7. Amount of times a sequence must be observed per sample to be retained for relative abundance calculation (recommended to set this >1 to remove singletons [highly likely to be PCR artifacts or sequencing errors])
8. Mimimum relative abundance as a percentage (eg. 2 or 0.1) for a sequence to be reported 
9. Number of threads to use for minimap2 mapping  

eg. `qontas.sh sample.fastq.gz ref.fa Qontas_Out sample 600 650 10 1 12`  

`qontas` will print these values to screen and wait 10 seconds before running to give you a chance to cancel if anything is wrong  

## TO DO
- [ ] automate pipeline in (better) single executable
- [ ] create test dataset  
- [ ] give option to retain or detele temp files
- [ ] modify to accept a list of input FASTQ files (TSV format) and output a single merged feature table  
   * alternatively, write a script that combines all the individual outputs 
   * Will need to modify to generate md5 read names so that they are groupable between samples  
- [ ] include flag to modify minimap -x flag allowing PacBio (`-x map-pb`) or Illumina (`-x sr`) reads (these will probably need to be merged beforehand)