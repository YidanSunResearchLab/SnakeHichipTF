import os
import argparse
import subprocess
import json
import re
import glob
import sys

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Generate genome indices and restriction enzyme fragments for the HiChIP pipeline.")
    parser.add_argument("-o", "--output_dir", required=True, help="Absolute Path to the output directory for storing generated genome files (default: ./genome/).")
    parser.add_argument("--organism", required=True, help="Reference genome assembly to generate (e.g., 'hg38' or 'mm10').")
    parser.add_argument("--organism_fasta", required=True, default="NA", help="Path to the organism FASTA file with 'chr' prefix (e.g., /path/to/hg38.fa). Default: NA.")
    parser.add_argument("--restriction_enzyme", default="mboi", help="Restriction enzyme to process. Options include 'mboi' (default), 'hindiii', 'dpnii', 'bglii', 'ncoi', 'msei', 'hinfI', 'mnase', 'arima'.")
    parser.add_argument("--bin_size", type=int, default=5000, help="Bin size in base pairs for fragment generation (default: 5000 bp).")
    parser.add_argument("--threads", type=int, default=30, help="Number of CPU threads to use for parallel processing (default: 30).")

    # Parse known arguments, capturing extra Snakemake args
    args, extra_args = parser.parse_known_args()

    # Ensure output directory exists
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Determine script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Check if the FASTA file exists
    organism = args.organism
    organism_fasta = args.organism_fasta
    fasta_path = f"{output_dir}/organisms/{organism}/{organism}.txt"
    if not os.path.exists(fasta_path):
        # If organism_fasta is not provided, raise an error
        if not os.path.exists(organism_fasta):
            sys.stderr.write(
                f"Dear user: The file '{fasta_path}' does not exist. "
                "To set up the genome index and restriction enzyme sites for "
                f"'{organism}', please specify --organism_fasta with the path to the organism's FASTA file.\n"
            )
            sys.exit(1)
        # If organism_fasta is provided, the pipeline will proceed (e.g., copy_fasta rule will handle it)

    # Build Snakemake command
    snakemake_cmd = [
        "snakemake",
        "--snakefile", os.path.join(script_dir, "../rules/Genomesetup.Snakefile"),
        "--printshellcmds",
        "--directory", output_dir,
        "--use-conda",
        "--conda-prefix", os.path.join(output_dir, ".snakemake_conda"),
        #"--use-singularity",
        #"--singularity-args", f"-B {output_dir} -B {singularity_tmpdir} -B {singularity_cachedir}",
        "--config",
        f"output_dir={args.output_dir}",
        f"organism={args.organism}",
        f"organism_fasta={args.organism_fasta}",
        f"restriction_enzyme={args.restriction_enzyme}",
        f"bin_size={args.bin_size}",
        f"threads={args.threads}",
        "--cores", str(args.threads),
        ]

    # Append any extra Snakemake arguments
    if extra_args:
        snakemake_cmd.extend(extra_args)

    # Generate unique log file name
    log_pattern = os.path.join(output_dir, f"genomesetup_run_{organism}_*.log")
    log_number = max([int(re.search(r"genomesetup_run_(.+)(\d+)\.log", f).group(2)) 
                    for f in glob.glob(log_pattern)], default=0) + 1
    log_file_path = os.path.join(output_dir, f"genomesetup_run_{organism}_{log_number}.log")

    # Log the command line
    command_line = " ".join(sys.argv)

    # Execute Snakemake with logging
    with open(log_file_path, "w") as log_file:
        # Write the command line to the log file
        log_file.write(f"Command line: {command_line}\n")
        log_file.write("-" * 80 + "\n")  # Separator for clarity
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