import subprocess
import argparse
import os.path
import pysam
from io import StringIO
import itertools
import numpy as np
import random
import re
from tempfile import TemporaryFile
import time
import copy
import logger
import sys
from shutil import copyfile

sam_columns = {"query_name": 0, "flag": 1, "chr_name": 2, "pos": 3, "mapq": 4,
               "cigar": 5, "mate_name": 6, "mate_pos": 7, "tlen": 8, "seq": 9, "quality": 10}

flag_table_proper = {(0, 16): (99, 147), (2048, 16): (99, 147), (0, 2064): (99, 147), (2048, 2064): (99, 147),  # + -
                     (16, 0): (83, 163), (16, 2048): (83, 163), (2064, 0): (83, 163), (2064, 2048): (83, 163),  # - +
                     (0, 0): (65, 129), (0, 2048): (65, 129), (2048, 0): (65, 129), (2048, 2048): (65, 129),    # + +
                     (16, 16): (113, 177), (16, 2064): (113, 177), (2064, 16): (113, 177), (2064, 2064): (113, 177)}  # - -

def filter_main(hicpro_bam, mapq, outdir, prefix, threads, optical_duplicate_distance, to_file=False):
    sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
    print(time.ctime() + " starting filtering operation with HiC-Pro BAM input")

    # Check if the input BAM file exists
    if not os.path.exists(hicpro_bam):
        exit(f"Error: Input BAM file {hicpro_bam} does not exist. Exiting!")
    if not hicpro_bam.endswith(".bam"):
        exit("Error: Input file must be a BAM file. Exiting!")
    if int(threads) < 1:
        exit("Error: Number of threads (-t) should be a positive integer. Exiting!")

    # Set filenames for intermediate and output files
    paired_filename, chromap_sorted_filename, combined_chromap_filename, qc_filename = set_filenames(outdir, prefix)

    # Open the HiC-Pro BAM file
    chromap = pysam.AlignmentFile(hicpro_bam, "rb")

    # Check if the BAM is sorted by query name; if not, sort it
    if not is_sorted_queryname(chromap.header):
        print(time.ctime() + " calling samtools sort for " + hicpro_bam + " storing in " + chromap_sorted_filename)
        pysam.sort("-o", chromap_sorted_filename, "-n", "-@", str(threads), hicpro_bam)
    else:
        copyfile(hicpro_bam, chromap_sorted_filename)

    # Merge (in this case, just use the sorted BAM since we have one input)
    print(time.ctime() + " using sorted HiC-Pro BAM as combined input")
    pysam.merge("-n", "-f", combined_chromap_filename, chromap_sorted_filename)

    # Fix mate information
    print(time.ctime() + " running samtools fixmate")
    pysam.fixmate("-m", "-@", str(threads), combined_chromap_filename, paired_filename + ".fixmated.bam")

    # Sort by coordinates
    print(time.ctime() + " sorting by coordinates")
    pysam.sort("-o", paired_filename + ".srt.bam", "-@", str(threads), paired_filename + ".fixmated.bam")

    # Mark and remove duplicates
    print(time.ctime() + " calling samtools markdup")
    try:
        proc = subprocess.run(
            " ".join([
                "samtools", "markdup",
                "-r", "-m", "s", "-s", "-f", paired_filename + ".fixmated.markdup.stats",
                "-@", str(threads), "-d", str(optical_duplicate_distance),
                paired_filename + ".srt.bam", paired_filename + ".rmdup.bam"
            ]),
            shell=True,
            check=True
        )
    except subprocess.CalledProcessError:
        print(time.ctime() + " modern samtools markdup failed, falling back to basic version")
        proc = subprocess.run(
            f"samtools markdup -r -s -@ {threads} {paired_filename}.srt.bam {paired_filename}.rmdup.bam 2> {paired_filename}.fixmated.markdup.stats",
            shell=True,
            check=True
        )

    # Filter for properly paired reads
    print(time.ctime() + " filtering unpaired reads after markdup")
    paired_only_bam = paired_filename + ".rmdup.paired.bam"
    pysam.view("-b", "-f", "2", "-o", paired_only_bam, paired_filename + ".rmdup.bam", catch_stdout=False)

    # Generate flagstat for QC
    print(time.ctime() + " calling samtools flagstat on mapped and duplicate-removed file")
    proc = subprocess.Popen("samtools flagstat " + paired_only_bam + " > " + paired_filename + ".rmdup.flagstat", shell=True)
    proc.communicate()

    # Final sort by query name
    print(time.ctime() + " sorting by query name for final output")
    pysam.sort("-o", paired_filename + ".srtn.rmdup.bam", "-@", str(threads), "-n", paired_only_bam)

    print(time.ctime() + " filtering complete")
    return paired_filename + ".srtn.rmdup.bam"

def set_filenames(outdir, prefix):
    tempdir = outdir + "/tempfiles"
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)
    paired_filename = outdir + "/" + prefix + ".paired"
    chromap_sorted_filename = tempdir + "/" + prefix + ".sorted.bam"
    combined_chromap_filename = tempdir + "/" + prefix + ".merged.srtn.bam"
    qc_filename = outdir + "/" + prefix + ".qc"
    return paired_filename, chromap_sorted_filename, combined_chromap_filename, qc_filename

def is_sorted_queryname(header):
    if "HD" in header:
        if "SO" in header["HD"]:
            if header["HD"]["SO"] == "queryname":
                return True
    return False
