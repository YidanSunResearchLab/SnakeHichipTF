#!/usr/bin/env python3
import os
import argparse
import subprocess
import json
import re
import glob
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Run TF scPrinter pipeline (supports BAM or fragments TSV(.gz) input) via Snakemake."
    )

    # Required core inputs (keep)
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Path to the output directory where results will be saved."
    )
    parser.add_argument(
        "--organism", required=True,
        help="Reference genome assembly (e.g., 'hg38' or 'mm10')."
    )
    parser.add_argument(
        "--genome_dir", default="",
        help="Absolute path to the directory containing pre-generated genome files."
    )

    # Keep/require enhancer bed input (explicit)
    parser.add_argument(
        "--enhancer_bedpe", required=True,
        help="Space-separated list of enhancer BEDPE files."
    )

    # require exactly one of BAM list or TSV(.gz) list
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--atacseq_bam",
        help="Space-separated list of ATAC-seq BAM files (e.g., 's1.bam s2.bam')."
    )
    group.add_argument(
        "--fragments_tsv",
        help="Space-separated list of fragments TSV.GZ files (e.g., 's1.tsv.gz s2.tsv.gz')."
    )

    # Samples must be explicit (recommended for reproducibility)
    parser.add_argument(
        "--samples", required=True,
        help="Space-separated list of sample names corresponding to the input files."
    )

    # Keep threads and the early required variables
    parser.add_argument(
        "--threads", type=int, default=30,
        help="Number of CPU threads to use (default: 30)."
    )

    # Mode (bulk vs single-cell) is important for downstream steps, so keep it explicit
    parser.add_argument(
        "--mode", default="bulk", choices=["bulk", "singlecell"],
        help="Analysis mode: 'bulk' for bulk ATAC-seq, 'singlecell' for single-cell ATAC-seq (default: 'bulk')."
    )


    # Parse known arguments, capturing extra Snakemake args
    args, extra_args = parser.parse_known_args()

    # Assign variables
    output_dir = os.path.abspath(args.output_dir)
    organism = args.organism
    genome_dir = os.path.abspath(args.genome_dir)
    enhancer_bedpe = args.enhancer_bedpe.split()
    samples = args.samples.split()
    mode = args.mode
    threads = args.threads

    # Input files: either bam list or fragments tsv list
    atacseq_bam = args.atacseq_bam.split() if args.atacseq_bam else []
    fragments_tsv = args.fragments_tsv.split() if args.fragments_tsv else []
    input_type = "bam" if atacseq_bam else "fragments_tsv"
    input_files = atacseq_bam if atacseq_bam else fragments_tsv

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Determine script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Validate inputs length
    if len(input_files) != len(samples):
        sys.stderr.write(
            f"Error: Number of input files ({len(input_files)}) does not match number of samples ({len(samples)}).\n"
            "Please ensure --samples matches the number of files you provide.\n"
        )
        sys.exit(1)

    if len(enhancer_bedpe) not in (1, len(samples)):
        sys.stderr.write(
            f"Error: --enhancer_bedpe must have either 1 file (shared) or the same number as samples ({len(samples)}).\n"
            f"Got {len(enhancer_bedpe)} enhancer_bedpe files.\n"
        )
        sys.exit(1)

    # Check if input files exist
    for f in input_files:
        if not os.path.exists(f):
            sys.stderr.write(f"Error: Input file '{f}' does not exist.\n")
            sys.exit(1)

    for bedpe in enhancer_bedpe:
        if not os.path.exists(bedpe):
            sys.stderr.write(f"Error: BEDPE file '{bedpe}' does not exist.\n")
            sys.exit(1)

    # Build a samplesheet-like mapping (to pass through config cleanly)
    # Each sample maps to its input file and (optionally) enhancer bedpe
    sample_map = {}
    for i, s in enumerate(samples):
        sample_map[s] = {
            "input": input_files[i],
            "enhancer_bedpe": enhancer_bedpe[i] if len(enhancer_bedpe) == len(samples) else enhancer_bedpe[0],
        }

    # Build Snakemake command
    snakemake_cmd = [
        "snakemake",
        "--snakefile", os.path.join(script_dir, "../rules/TFscprinter.Snakefile"),
        "--printshellcmds",
        "--directory", output_dir,
        "--use-conda",
        "--conda-prefix", os.path.join(genome_dir, ".snakemake_conda"),
        "--config",
        f"organism={organism}",
        f"genome_dir={genome_dir}",
        f"enhancer_bedpe={enhancer_bedpe}",
        f"input_type={input_type}",                 # "bam" or "fragments_tsv"
        f"samples={json.dumps(samples)}",           # keep for backward compatibility
        f"sample_map={json.dumps(sample_map)}",     # preferred structured mapping
        f"mode={mode}",
        f"threads={threads}",
        "--cores", str(threads),
    ]

    # Append any extra Snakemake arguments
    if extra_args:
        snakemake_cmd.extend(extra_args)

    # Generate unique log file name
    log_pattern = os.path.join(output_dir, "tfscprinter_run_*.log")
    log_number = max(
        [int(re.search(r"tfscprinter_run_(\d+)\.log", f).group(1)) for f in glob.glob(log_pattern)],
        default=0
    ) + 1
    log_file_path = os.path.join(output_dir, f"tfscprinter_run_{log_number}.log")

    # Log the command line
    command_line = " ".join(sys.argv)

    # Execute Snakemake with logging
    with open(log_file_path, "w") as log_file:
        log_file.write(f"Command line: {command_line}\n")
        log_file.write("-" * 80 + "\n")
        process = subprocess.Popen(
            snakemake_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in process.stdout:
            print(line, end="")
            log_file.write(line)
        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, snakemake_cmd)


if __name__ == "__main__":
    main()