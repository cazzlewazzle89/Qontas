# Qontas

### ✈️ Quantification of ONT Amplicon Sequences ✈️

**QONTAS** is a pipeline designed to bin, clip, and denoise Nanopore amplicon data. By using a "Cookie Cutter" approach, it extracts specific genomic regions defined in a BED file, flips them to a consistent orientation, and uses VSEARCH to collapse noisy reads into high-confidence biological variants.

## Key Features

* **Hard Clipping:** Strictly extracts sequences within BED coordinates to remove primer/adapter/off-target noise.
* **Orientation Recovery:** Automatically reorients reverse-mapped reads to the forward strand.
* **Denoising:** Uses `vsearch` dereplication to filter out random sequencing errors (singletons).
* **Storage Optimized:** Automatically cleans up large intermediate BAM files after processing.

## SETUP

### 1. Clone the Repository

```bash
git clone https://github.com/cazzlewazzle89/Qontas.git ~/Software/QONTAS
chmod +x ~/Software/QONTAS/qontas.py
```

### 2. Create Conda Environment

We recommend using `mamba` for faster dependency resolution:
```bash
mamba create -n qontas -c bioconda -c conda-forge minimap2 samtools vsearch pysam pandas biopython
conda activate qontas
```

## USAGE

The pipeline is run via `qontas.py`. It requires a BED file that includes your target length constraints in columns 5 and 6.

### Amplicon Panels (Multi-Region Support)

QONTAS can also handle multi-target amplicon panels in a single run. You can provide as many regions as you like in the BED file. The script will:
1.  Bin reads to each region simultaneously using a global alignment.
2. Process each region independently.
3.  Generate separate FASTA and TSV outputs for every region listed.

### BED File Format

Ensure your BED file is tab-separated and follows this structure:

```text
#Chrom   Start   End     Name    MinLen   MaxLen
chr1     1658    1717    rpoB    600      650
chr1     2500    3100    gyrA    580      620
```

### Running the Script

```bash
python qontas.py \
    -i sample.fastq.gz \
    -r reference.fasta \
    -b regions.bed \
    -s Sample_Name \
    -o Qontas_Out \
    --mincount 10
```

### Arguments:

| Flag | Description |
| --- | --- |
| `-i` | Input FASTQ (can be gzipped). |
| `-r` | Reference FASTA for the whole genome or target regions. |
| `-b` | BED file defining regions, names, and length filters. |
| `-s` | **Sample Name**: Used to create a unique subdirectory. |
| `-o` | **Output Directory**: Parent folder for all results. |
| `--mincount` | Minimum identical reads required to retain a cluster (Default: 10). |

## OUTPUT STRUCTURE

QONTAS creates a clean, sample-centric directory structure:

```text
Qontas_Out/
└── Sample_Name/
    ├── rpoB.fa    # Multi-FASTA of denoised variants for rpoB
    ├── rpoB.tsv   # TSV with counts/abundance for rpoB
    ├── gyrA.fa    # Multi-FASTA of denoised variants for gyrA
    └── gyrA.tsv   # TSV with counts/abundance for gyrA
```

## TO DO

* [ ] **Variant Calling:** Add a step to generate SNV strings (e.g., A26G) relative to reference.
* [ ] **Multi-Sample:** Add a script to merge all `.tsv` files into a single feature table.
* [ ] **Test Data:** Provide a small subset of reads for benchmarking.
