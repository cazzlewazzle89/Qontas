# Qontas

### ✈️ Quantification of ONT Amplicon Sequences ✈️

**QONTAS** (Quantification of ONT Amplicon Sequences) is a tool designed to take Nanopore amplicon sequencing reads and generate a clean table of unique variants, their counts, and relative abundances. It uses a "Cookie Cutter" approach to bin, clip, and denoise data based on specific genomic regions of interest.

## Key Features

* **Per-Region Length Filtering:** Filters reads based on their original length before clipping, allowing for diverse amplicon panels.
* **Hard Clipping:** Extracts exact genomic windows to remove primer and adapter artifacts.
* **Strand-Agnostic Clustering:** Automatically merges forward and reverse-complement reads into single biological clusters.
* **Storage Optimized:** Automatically removes large intermediate BAM files after processing to save disk space.

## SETUP

### 1. Clone the Repository

```bash
git clone https://github.com/cazzlewazzle89/Qontas.git ~/Software/QONTAS
chmod +x ~/Software/QONTAS/qontas.py
```

### 2. Create Conda Environment

```bash
mamba create -n qontas -c bioconda -c conda-forge minimap2 samtools vsearch pysam pandas biopython
conda activate qontas
```

## USAGE

QONTAS is run using the `qontas.py` script. It supports multi-region amplicon panels via a single BED file.

### BED File Format

The BED file should be tab-separated. Columns 5 and 6 are used to define the expected length of the **original read** for that specific target. For example:

```text
#Chrom   Start   End     Name    MinLen   MaxLen
chr1     1658    1717    rpoB    600      650
chr1     2500    3100    gyrA    1100     1250
```

### Running the Script

```bash
python qontas.py -i sample.fastq.gz -r ref.fa -b regions.bed -s SampleName -o Qontas_Out
```

### Command Line Arguments:

| Flag | Name | Description | Default |
| --- | --- | --- | --- |
| `-i` | Input | Sequencing reads in FASTQ format (can be gzipped) | Required |
| `-r` | Reference | Reference sequence in FASTA format | Required |
| `-b` | BED | BED file defining regions and length constraints | Required |
| `-s` | Sample | Sample name used for subdirectory and file naming | Required |
| `-o` | Outdir | Parent directory for results | `Qontas_Out` |
| `-c` | Mincount | Minimum reads per cluster to be retained | `10` |
| `-t` | Threads | Number of threads for Minimap2 mapping | `4` |

## OUTPUT STRUCTURE

For each sample, QONTAS generates a dedicated folder containing two files per region defined in the BED file:

```text
Qontas_Out/
└── SampleName/
    ├── rpoB.fa    # Multi-FASTA of denoised variants
    └── rpoB.tsv   # TSV with ClusterID, Count, and Relative Abundance
    ├── gyrA.fa    # Multi-FASTA of denoised variants
    └── gyrA.tsv   # TSV with ClusterID, Count, and Relative Abundance
```

## TO DO

* [ ] **Variant Calling:** implement SNV string generation (e.g., S531L) relative to reference.
* [ ] **Multi-Sample:** Script to merge individual `.tsv` files into a single master feature table.
