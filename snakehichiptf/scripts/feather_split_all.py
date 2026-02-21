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
downsample_pairs = 0
chip_peaks = None
optical = 0

def main():
    parser = setup_arg_parsers()
    parse_args(parser)
    sys.stdout = logger.Logger(outdir + "/" + prefix + ".feather.log")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    from feather_split_rongxin import split_main
    split_main(split_filename, outdir, prefix, length_threshold, per_chr_bedpe, generate_hic, downsample_pairs, chip_peaks)

def parse_args(parser):
    args = parser.parse_args()
    for arg in args.__dict__:
        val = args.__dict__[arg]
        if val is not None:
            globals().update({arg: val})

def setup_arg_parsers():
    parser = argparse.ArgumentParser(description="Split a BAM file by contact length threshold into BEDPE files.")
    parser.add_argument("-o", "--outdir", help="Path to the output directory.", required=True)
    parser.add_argument("-p", "--prefix", help="Prefix for output files.")
    set_split_parser_args(parser)
    return parser

def set_split_parser_args(parser):
    parser.add_argument("-l", "--length-threshold", help="Length threshold for long/short splitting.", type=int)
    parser.add_argument("-s", "--split-filename", help="BAM or SAM file to split.")
    parser.add_argument("-c", "--per-chr-bedpe", help="Generate one BEDPE file per chromosome.", type=str2bool, nargs='?', const=True, default="t")
    parser.add_argument("-j", "--generate-hic", help="Generate a .hic file.", type=str2bool, nargs='?', const=True, default=generate_hic)
    parser.add_argument("-d", "--downsample_pairs", help="Downsample size", default="0")
    parser.add_argument("-a", "--chip-peaks", help="1D anchors BED file (e.g., ChIP/ATAC-seq peaks).", default=chip_peaks)

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

main()