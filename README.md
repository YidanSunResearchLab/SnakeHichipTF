# SnakeHichipTF

**A Comprehensive Snakemake Pipeline for integrating HiChIP / PLAC-Seq Processing and Transcription Factor Binding Analysis**

---

## Overview

**SnakeHichipTF** is a modular, reproducible, and scalable Snakemake-based workflow for integrative analysis of **HiChIP** and **PLAC-seq** datasets, with optional integration of **ATAC-seq–based transcription factor (TF) footprinting**.

![Image](docs/overview.png)


The pipeline supports:

* End-to-end processing from FASTQ → interaction maps → loop calling → TF binding inference
* Multiple mapping engines: **HiC-Pro**, **Chromap**
* Multiple loop callers: **MAPS**, **FitHiChIP**, **HiChipper**, **HiCDCPlus**
* AI-based TF binding inference via **scPrinter**
* Footprinting-based TF analysis using **TOBIAS**
* Genome setup utilities (restriction fragments, indices)
* Read downsampling and QC reporting
* Fully reproducible execution via **Snakemake**

SnakeHichipTF is designed for flexible experimental designs, including bulk HiChIP, PLAC-seq, and integration with matched ATAC-seq data.

---

## Key Features

### 1. Modular Architecture

Each analytical stage is implemented as an independent Snakemake rule:

* Genome preparation
* Significant interaction analysis
* Differential interaction analysis
* TF binding and footprint analysis

This allows selective execution and easy customization.

---

### 2. Multi-Algorithm Loop Calling for HiChIP/PLAC-seq 

Supported loop callers:

* **MAPS**
* **FitHiChIP**
* **HiChipper**
* **HiCDCPlus**

Users can benchmark across algorithms within a single workflow.

---

### 3. AI-Based TF Binding Inference (scPrinter)

SnakeHichipTF integrates **scPrinter**, enabling:

* seq2PRINT model training
* TF binding signal inference
* Motif scanning and TFBS scoring
* Comparative TF activity visualization

---

### 4. ATAC-Seq Footprinting (TOBIAS)

For chromatin accessibility-based TF analysis:

* Bias correction
* Footprint detection
* Motif enrichment
* TF activity scoring

---

### 5. Reproducibility and Scalability

* Snakemake-managed dependency control
* Conda environment integration
* HPC cluster compatibility
* Parallelized execution
* Resume capability

---

# Installation

## 1. Clone Repository

```bash
git clone https://github.com/YidanSunResearchLab/SnakeHichipTF.git
cd SnakeHichipTF
```

## 2. Create Environment

```bash
conda create -n SnakeHichipTF -c conda-forge python>=3.8
conda activate SnakeHichipTF
pip install .
```

(Advanced users may use a dedicated environment for GPU-enabled scPrinter workflows.)

## 3. Verify Installation

```bash
Snakehichip --help
```

---

# Quick Start Guide

SnakeHichipTF provides several high-level commands.

---

## 1️⃣ Genome Setup

Prepare genome indices and restriction fragments.

```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /path/to/genome.fa \
  -o ./genome \
  --restriction_enzyme mboi \
  --bin_size 5000 \
  --threads 32
```

### Outputs:

* Indexed genome
* Restriction fragment BED
* Precomputed binning files

---

## 2️⃣ HiChIP / PLAC-Seq Processing

Process FASTQ files into interaction matrices and loop calls.

```bash
Snakehichip \
  -i /path/to/input_fastqs/ \
  --genome_dir /path/to/genome \
  -o /path/to/output_directory/ \
  --organism hg38 \
  --bin_size 5000 \
  --threads 32 \
  --downsample_size 50000000
```

### Steps Performed:

1. Read alignment
2. Filtering & deduplication
3. Contact matrix construction
4. Loop calling (multiple algorithms)
5. Quality control reporting

---

## 3️⃣ Differential HiChIP Analysis

```bash
Diffhichip \
  -s /path/to/sample_sheet.tsv \
  -o /path/to/output_directory/ \
  --threads 16 \
  --bin_size 5000 \
  --fitType mean \
  --fdr 0.05
```

### Inputs:

* Sample sheet describing sample file paths and condition labels, examples are in the data/ folder

### Outputs:

* Differential interaction results
* FDR-filtered loop sets

---

## 4️⃣ TF Binding Analysis via scPrinter

Run AI-based TF binding inference using ATAC-seq input and loop regions.

```bash
TFscprinter \
  -o /path/to/output_directory/ \
  --organism hg38 \
  --atacseq_bam "/path/to/sample.bam" \
  --enhancer_bedpe "/path/to/significant_loops.bedpe" \
  --samples "Sample1" \
  --threads 32
```

### Outputs:

* seq2PRINT trained models
* TFBS bigWig tracks
* Motif scanning results
* TFBS score tables
* Scatter plot comparisons

---

## 5️⃣ ATAC-Seq Footprinting with TOBIAS

```bash
TFtobias \
  -o /path/to/output_directory/ \
  --organism hg38 \
  --atacseq_bam "/path/to/sample1.bam /path/to/sample2.bam" \
  --enhancer_bedpe "/path/to/sample1.bedpe /path/to/sample2.bedpe" \
  --samples "Sample1 Sample2" \
  --threads 32 \
  --motif "/path/to/JASPAR2024_CORE_non-redundant_pfms.meme"
```

### Outputs:

* Corrected bigWig tracks
* Footprint regions
* Motif enrichment reports
* TF activity scores

---

# Directory Structure

```
SnakeHichipTF/
│
├── snakehichiptf/
│   ├── workflows/
│   ├── rules/
│   ├── scripts/
│
├── docs/
│   ├── Genomesetup.md
│   ├── Snakehichip.md
│   ├── TFscprinter.md
│   ├── TFtobias.md
│
└── setup.py
```

---

# Documentation

Detailed documentation for each module:

- [Genomesetup.py](docs/Genomesetup.md)
- [Snakehichip.py](docs/Snakehichip.md)
- [TFscprinter.py](docs/TFscprinter.md)
- [TFtobias.py](docs/TFtobias.md)

See the `docs/` directory for command-level explanations.

---

# Recommended Usage Strategy

For large-scale projects:

1. Perform genome setup once.
2. Run Snakehichip for each sample.
3. Aggregate and run Diffhichip if comparing conditions.
4. Perform TFscprinter or TFtobias for mechanistic interpretation.

---

# Contributing

Issues and pull requests are welcome:

[https://github.com/YidanSunResearchLab/SnakeHichipTF](https://github.com/YidanSunResearchLab/SnakeHichipTF)

---

# License

See [LICENSE](LICENSE) for details.
