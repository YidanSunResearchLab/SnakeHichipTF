# HiChIPSnake: A Comprehensive Pipeline for HiChIP/PLAC-Seq Data Analysis

**HiChIPSnake** is a Snakemake-based pipeline designed for end-to-end analysis of HiChIP and PLAC-seq datasets. It streamlines the processing of chromosome conformation data by integrating robust mapping tools such as **HiC-Pro** and **Chromap**, and supports multiple state-of-the-art loop-calling algorithms including **MAPS**, **FitHiChIP**, **HiChipper**, and **HiCDCPlus**. Additionally, **TOBIAS** is incorporated for transcription factor footprinting analysis.

HiChIPSnake also includes features for **genome setup**, **read downsampling**, and **comprehensive quality control**, offering researchers a reproducible, modular, and scalable solution for analyzing chromatin interaction data across diverse experimental designs.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/YidanSunResearchLab/HichipSnake.git
   cd HichipSnake
   ```

2. **Install Dependencies**:
   - Set up the environment:
  ```bash
  conda create -n HichipSnake -c conda-forge python>=3.8
  conda activate HichipSnake
  pip install .
  ```

3. **Verify Installation**:
   ```bash
   Hichipsnake --help
   ```

## Quick Start

HiChIPSnake offers three main commands:

### 1. Genome Setup
Prepare genome indices and restriction fragments:
```bash
Genomesetup \
  --organism hg38 \
  --organism_fasta /path/to/your/genome.fa \
  -o ./genome \
  --restriction_enzyme mboi \
  --bin_size 5000 \
  --threads 100
```

### 2. HiChIP/PLAC-Seq Processing
Process FASTQs into interaction maps:
```bash
Hichipsnake \
  -i /home/syidan/syidan/Projects/hichip_real_samples/ \
  --genome_dir /abs/path/to/genome \
  -o /home/syidan/syidan/Projects/hichip_real_samples_output/ \
  --threads 100 \
  --organism hg38 \
  --bin_size 5000 \
  --downsample_size 50000000
```

### 3. ATAC-Seq Footprinting
Analyze ATAC-seq data for footprints:
```bash
ATACfootprinting \
  -o /home/syidan/syidan/Projects/hichip_test_samples_output \
  --organism hg38 \
  --atacseq_bam "/path/to/sample1.bam /path/to/sample2.bam" \
  --enhancer_bed "/path/to/sample1.bedpe /path/to/sample2.bedpe" \
  --samples "Sample1 Sample2" \
  --threads 100 \
  --motif "/Path/to/JASPAR2024_CORE_non-redundant_pfms_meme03142025"
```

## Documentation
For detailed usage, see:
- [Genomesetup.py](docs/Genomesetup.md)
- [Hichipsnake.py](docs/Hichipsnake.md)
- [ATACfootprinting.py](docs/ATACfootprinting.md)

## Contributing
Feel free to open issues or submit pull requests on [GitHub](https://github.com/YidanSunResearchLab/HiChIPSnake).

## License
See [LICENSE](LICENSE) for details.


