# TFscprinter: AI-Based Transcription Factor Binding Inference

`TFscprinter` performs transcription factor (TF) binding inference using the **scPrinter (seq2PRINT)** framework.

It supports both:

* Bulk ATAC-seq (BAM input)
* Pre-generated fragments files (`.tsv.gz`)
* Single-cell ATAC-seq (optional mode)

Enhancer regions (typically derived from HiChIP/PLAC-seq loop calling) are used as candidate regulatory regions for TF activity analysis.

---

## Usage

```bash
TFscprinter [options]
```

---

## Required Arguments

* `-o/--output_dir`
  Output directory where results will be written.

* `--organism`
  Reference genome assembly (e.g., `hg38`, `mm10`). Must match genome used in `Genomesetup`.

* `--enhancer_bedpe`
  Space-separated list of enhancer BEDPE files.
  These are typically significant loop interactions from the HiChIP/PLAC-seq step.

* `--samples`
  Space-separated sample names corresponding to input files.

* **One of the following input types (mutually exclusive):**

  * `--atacseq_bam`
    Space-separated list of filtered ATAC-seq BAM files
    (e.g., `"sample1.filtered.bam sample2.filtered.bam"`)

  * `--fragments_tsv`
    Space-separated list of fragment files (`.tsv.gz`)
    (e.g., `"sample1.fragments.tsv.gz sample2.fragments.tsv.gz"`)

---

## Optional Arguments

* `--genome_dir`
  Path to genome directory generated during `Genomesetup`.

* `--threads`
  Number of CPU threads to use.
  Default: `30`

* `--mode`
  Analysis mode:

  * `bulk` (default)
  * `singlecell`

* `--plot_aggregate`
  Optional parameters passed to aggregate plotting.

---

## Input Validation Rules

* Number of input files must equal number of `--samples`
* `--enhancer_bedpe` must contain:

  * either 1 file (shared across samples), or
  * same number of files as samples
* All input files must exist
* BAM files must be:

  * filtered
  * coordinate-sorted
  * indexed

---

## Example (Bulk ATAC-seq Bam)

```bash
TFscprinter \
  -o /path/to/output_directory/ \
  --organism hg38 \
  --genome_dir /path/to/genome \
  --atacseq_bam "sample1.filtered.bam sample2.filtered.bam" \
  --enhancer_bedpe "loops.significantInteractions.bedpe" \
  --samples "Sample1 Sample2" \
  --threads 32 \
  --mode bulk
```

---

## Example (Fragments Input)

```bash
TFscprinter \
  -o /path/to/output_directory/ \
  --organism hg38 \
  --fragments_tsv "sample1.fragments.tsv.gz sample2.fragments.tsv.gz" \
  --enhancer_bedpe "loops.significantInteractions.bedpe" \
  --samples "Sample1 Sample2" \
  --threads 32
```

---


# TFscprinter Workflow Overview

For each sample, TFscprinter performs:

1. Fragment import and filtering
2. Peak calling (seq2PRINT preset)
3. Model configuration generation
4. seq2PRINT model training
5. TF binding score inference
6. Motif scanning
7. TFBS score extraction
8. Comparative visualization

---

# TFscprinter Output Structure

```
output_directory/
├── fragments/
│   └── *.fragments.tsv.gz
│
├── seq2print/
│   ├── ATAC_scprinter.h5ad
│   ├── seq2print_cleaned_narrowPeak.bed
│   ├── regions.bed
│   ├── regions_filt.bed
│   ├── configs/
│   │   └── ATAC_<sample>_fold0.JSON
│   ├── model/
│   │   └── ATAC_<sample>*.pt
│   ├── tfbs_<sample>/
│   │   └── dummy_TFBS.bigwig
│   ├── TFBS_scores.csv
│   ├── TFBS_scatter.png
│   └── scprinter_done.txt
│
└── tfscprinter_run_<n>.log
```

---

## Key Output Files

### `model/`

Trained seq2PRINT models for each sample.

### `TFBS_scores.csv`

Motif-level TF binding scores across enhancer regions.

### `TFBS_scatter.png`

Pairwise TF binding comparison plots between samples.

### `tfbs_<sample>/`

BigWig tracks containing inferred TF binding signal.

---

# Notes

* GPU is recommended for seq2PRINT training.
* Ensure CUDA environment is properly configured.
* For reproducibility, sample names must match input file order.
* Enhancer BEDPE files are typically derived from significant_interactions/ folder:

  * MAPS
  * FitHiChIP
  * HiChipper
  * HiCDCPlus

---

# Log Files

Each run generates a uniquely numbered log file:

```
tfscprinter_run_<n>.log
```

This file captures:

* Full command line
* Snakemake execution output
* Error messages (if any)

---
