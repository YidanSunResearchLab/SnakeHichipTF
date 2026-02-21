#!/usr/bin/env python

import argparse
import subprocess
import pysam
import pybedtools
import os
from tempfile import TemporaryFile
import sys
import logger

threads = 2
mapq = 30
fastq1 = ""
fastq2 = ""
chromap_index = ""
outdir = ""
genome = ""
re_pattern = ""
prefix = ""
length_threshold = 1000
split_filename = ""
per_chr_bedpe = "1"
generate_hic = "1"
chip_peaks = None
optical = 0

def main():
    parser = setup_arg_parsers()
    parse_args(parser)
    sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    from feather_filter_chr import filter_main
    filter_main(fastq1, fastq2, chromap_index, chromap_index_fa, mapq, outdir, prefix, threads, optical, to_file=True)

def parse_args(parser):
    args = parser.parse_args()
    for arg in args.__dict__:
        val = args.__dict__[arg]
        if val is not None:
            globals().update({arg: val})

def setup_arg_parsers():
    parser = argparse.ArgumentParser(description="Run chromap mapping and filter paired-end reads, generating a filtered BAM file.")
    parser.add_argument("-o", "--outdir", help="Path to the output directory.", required=True)
    parser.add_argument("-p", "--prefix", help="Prefix for output files.")
    set_filter_parser_args(parser)
    return parser

def set_filter_parser_args(parser):
    parser.add_argument("-f1", "--fastq1", help="First FASTQ or BAM file.", required=True)
    parser.add_argument("-f2", "--fastq2", help="Second FASTQ or BAM file.", required=True)
    parser.add_argument("-b", "--chromap-index", help="Genome index for chromap.")
    parser.add_argument("-fa", "--chromap-index-fa", help="Genome FASTA file for chromap.")
    parser.add_argument("-q", "--mapq", help="Minimum MAPQ threshold.", type=int)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int)
    parser.add_argument("-d", "--optical-duplicate-distance", help="Distance for optical duplicate removal.", type=int)

main()