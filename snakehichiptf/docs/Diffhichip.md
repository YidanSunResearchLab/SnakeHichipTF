# DiffHiChIP: Differential HiChIP Interaction Analysis

`DiffHiChIP` performs differential chromatin interaction analysis using within a Snakemake workflow.

This pipeline compares HiChIP interaction (BEDPE/BEDPE.GZ) files across biological conditions and identifies statistically significant differential loops.

---

## Usage

```bash
DiffHiChIP [options]
```

---

## Required Arguments

* `-s/--samples_sheet`
  Path to a **TSV samplesheet**.

  The samplesheet must contain two columns:

  ```
  samplename    condition
  ```

  * **Column 1:** Absolute path to per-sample all interaction txt.gz file
  * **Column 2:** Condition label (e.g., Control, Disease)

  Example:

  ```
  /path/sample1.txt.gz    Control
  /path/sample2.txt.gz    Control
  /path/sample3.txt.gz    Disease
  /path/sample4.txt.gz    Disease
  ```

* `-o/--output_dir`
  Output directory where Snakemake will run and results will be generated.

---

## Optional Arguments

* `--threads`
  Total cores passed to `snakemake --cores`.
  Default: `2`

* `--fdr`
  False discovery rate threshold for differential interactions.
  Default: `0.05`

* `--bin_size`
  Genomic bin size used in modeling.
  Default: `5000`

* `--fitType`
  Statistical fit type used in HiCDCPlus modeling.
  Default: `mean`

* `--conda_prefix`
  Custom conda environment prefix for Snakemake.

* Any additional Snakemake arguments can be passed after `--`, for example:

```bash
--rerun-incomplete --keep-going --jobs 200
```

---

## Example

```bash
DiffHiChIP \
  -s samplesheet.tsv \
  -o diffhichip_results \
  --threads 32 \
  --fdr 0.01 \
  --bin_size 10000 \
  --fitType mean \
  --conda_prefix /path/to/conda_envs \
  --rerun-incomplete --keep-going
```

---

## Workflow Overview

1. **Input validation**

   * Verifies samplesheet format
   * Confirms all listed BEDPE files exist

2. **Snakemake execution**

   * Uses `--use-conda` to ensure reproducible environments
   * Runs HiCDCPlus-based differential modeling
   * Generates differential interaction tables

3. **Statistical modeling**

   * Interaction counts are modeled across conditions
   * Dispersion and mean-variance trends are estimated
   * Differential loops are identified based on FDR threshold

---

## Output Structure

```
diffhichip_results/
├── differential_interactions.tsv
├── QC/
│   ├── dispersion_plot.pdf
│   ├── MA_plot.pdf
│   └── volcano_plot.pdf
├── diffhichip_run_1.log
└── .snakemake/
```

### Key Outputs

* **`differential_interactions.tsv`**
  Significant differential chromatin loops.

* **QC plots**
  Model diagnostics and statistical summaries.

---

## Logging

Each run generates an incremental log file:

```
diffhichip_run_1.log
diffhichip_run_2.log
...
```

The log captures:

* Full command line
* Snakemake execution output
* Errors and warnings

---

## Notes

* Input txt files must be:

  * Coordinate-sorted
  * Consistently formatted
  * Generated from the same HiChIP processing pipeline

* Replicates per condition are strongly recommended for reliable dispersion estimation.

* For cluster execution, additional Snakemake options (e.g., `--cluster`, `--jobs`) can be appended after `--`.

---
