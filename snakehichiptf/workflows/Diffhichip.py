#!/usr/bin/env python3
import os
import argparse
import subprocess
import glob
import re
import sys
from datetime import datetime

def read_samplesheet_paths(samplesheet: str):
    """
    Read first column (file path) from a TSV samplesheet with header:
    samplename    condition
    Return list of paths.
    """
    paths = []
    with open(samplesheet, "r") as f:
        header = next(f, None)
        if header is None:
            raise ValueError(f"Empty samplesheet: {samplesheet}")
        for line in f:
            line = line.strip()
            if not line:
                continue
            paths.append(line.split("\t")[0])
    return paths

def next_log_path(output_dir: str, prefix: str = "diffhichip_run_"):
    pattern = os.path.join(output_dir, f"{prefix}*.log")
    existing = glob.glob(pattern)
    max_n = 0
    for fp in existing:
        m = re.search(rf"{re.escape(prefix)}(\d+)\.log$", os.path.basename(fp))
        if m:
            max_n = max(max_n, int(m.group(1)))
    return os.path.join(output_dir, f"{prefix}{max_n + 1}.log")

def main():
    parser = argparse.ArgumentParser(
        description="Run DiffHiChIP Snakemake pipeline (HiCDCPlus differential interactions)."
    )
    parser.add_argument(
        "-s", "--samples_sheet", required=True,
        help="Path to samplesheet TSV. First column must be per-sample bedpe(.gz) paths; second column is condition."
    )
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Output directory where Snakemake will run."
    )
    parser.add_argument(
        "--snakefile", default=None,
        help="Path to Diffhichip.Snakefile. Default: <this_script_dir>/../rules/Diffhichip.Snakefile"
    )
    parser.add_argument("--threads", type=int, default=2, help="Total cores for snakemake --cores (default: 2).")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold passed to pipeline (default: 0.05).")
    parser.add_argument("--bin_size", type=int, default=5000, help="Bin size passed to pipeline (default: 5000).")
    parser.add_argument("--fitType", default="mean", help="fitType passed to pipeline (default: mean).")
    parser.add_argument(                                              # <-- 加这三行
        "--consensus_bedpe", default="",
        help="(optional) consensus BEDPE from consensus_interactions.py"
    )
    parser.add_argument(
        "--conda_prefix", default=None,
        help="Conda prefix for snakemake environments. If provided, passed via --conda-prefix."
    )
    parser.add_argument(
        "--dryrun", action="store_true",
        help="Run snakemake with -n (dry-run)."
    )

    # Capture any extra args after '--' and forward to snakemake
    args, extra_args = parser.parse_known_args()

    samplesheet = os.path.abspath(args.samples_sheet)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(samplesheet):
        sys.stderr.write(f"ERROR: samples_sheet not found: {samplesheet}\n")
        sys.exit(1)

    # Determine snakefile
    script_dir = os.path.dirname(os.path.abspath(__file__))
    snakefile = args.snakefile or os.path.join(script_dir, "../rules/Diffhichip.Snakefile")
    snakefile = os.path.abspath(snakefile)
    if not os.path.exists(snakefile):
        sys.stderr.write(f"ERROR: snakefile not found: {snakefile}\n")
        sys.exit(1)

    # Pre-parse samplesheet and check files exist
    try:
        sample_paths = read_samplesheet_paths(samplesheet)
    except Exception as e:
        sys.stderr.write(f"ERROR: failed to parse samplesheet: {e}\n")
        sys.exit(1)

    if len(sample_paths) == 0:
        sys.stderr.write("ERROR: samplesheet contains no sample paths.\n")
        sys.exit(1)

    missing = [p for p in sample_paths if not os.path.exists(p)]
    if missing:
        sys.stderr.write("ERROR: The following sample files listed in samplesheet do not exist:\n")
        for p in missing:
            sys.stderr.write(f"  - {p}\n")
        sys.exit(1)

    print("Detected sample files from samplesheet (first column):")
    for p in sample_paths:
        print(f"  {p}")

    # Build snakemake command
    snakemake_cmd = [
        "snakemake",
        "--snakefile", snakefile,
        "--directory", output_dir,
        "--printshellcmds",
        "--use-conda",
        "--cores", str(args.threads),
        "--config",
        f"samples_sheet={samplesheet}",
        f"output_dir={output_dir}",
        f"threads={args.threads}",
        f"fdr={args.fdr}",
        f"bin_size={args.bin_size}",
        f"fitType={args.fitType}",
        f"consensus_bedpe={args.consensus_bedpe}",
    ]

    if args.conda_prefix:
        snakemake_cmd.extend(["--conda-prefix", os.path.abspath(args.conda_prefix)])

    if args.dryrun:
        snakemake_cmd.append("-n")

    # forward extra snakemake args (e.g. --rerun-incomplete, --keep-going, --jobs 200, --cluster ...)
    if extra_args:
        snakemake_cmd.extend(extra_args)

    # Logging
    log_path = next_log_path(output_dir)
    command_line = " ".join(sys.argv)
    start_ts = datetime.now().isoformat(timespec="seconds")

    print("\nRunning command:")
    print(" ".join(snakemake_cmd))
    print(f"\nLog: {log_path}\n")

    with open(log_path, "w") as logf:
        logf.write(f"Start: {start_ts}\n")
        logf.write(f"Command line: {command_line}\n")
        logf.write("-" * 80 + "\n")
        process = subprocess.Popen(
            snakemake_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in process.stdout:
            print(line, end="")
            logf.write(line)
        process.wait()

    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, snakemake_cmd)

if __name__ == "__main__":
    main()
