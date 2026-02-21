#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
# Author: Yunjiang Qiu <serein927@gmail.com>
# Modified by: Armen Abnousi and optimization by ChatGPT
# File: restriction_cut.py
# Create Date: 2015-02-27 16:46:15
#########################################

import sys
import argparse
from Bio import SeqIO
import re
try:
    # Try to import the new function
    from Bio.SeqUtils import gc_fraction
except ImportError:
    # Fallback for older versions
    from Bio.SeqUtils import GC
    def gc_fraction(seq):
        return GC(seq) / 100.0


def find_site(fasta, seq, outfile, pos, binsize):
    with open(outfile, 'w') as f:
        # Determine if multiple enzyme recognition sequences are provided
        if (len(seq.split(",")) > 1):
            seqs = seq.split(",")
            sizes = [len(s) for s in seqs]
            poses = pos.split(",")
            poses = [int(x) for x in poses]
            mnase = False
            # Precompile patterns for each recognition sequence
            patterns = [re.compile(s.lower()) for s in seqs]
        else:
            # For a single recognition sequence, check if it's "mnase"
            mnase = True if seq.lower() == "mnase" else False
            seqs = [seq]
            sizes = [len(seq)]
            poses = [int(pos)] if not mnase else [0]
            patterns = [re.compile(seq.lower())] if not mnase else None

        # Process each sequence record in the FASTA file
        for seq_record in SeqIO.parse(fasta, "fasta"):
            if not mnase:
                sys.stderr.write("processing " + seq_record.id + "\n")
                record_seq = str(seq_record.seq).lower()  # Convert once per record
                site = []
                site_pos = []
                site_size = []
                for i, pattern in enumerate(patterns):
                    # Find all occurrences of the pattern in record_seq
                    current_site = [m.start() + 1 for m in pattern.finditer(record_seq)]
                    site.extend(current_site)
                    site_pos.extend([poses[i]] * len(current_site))
                    site_size.extend([sizes[i]] * len(current_site))
                if len(site) == 0:
                    continue
                # Sort the sites and associated data
                sorted_sites = sorted(zip(site, site_pos, site_size), key=lambda x: x[0])
                site, site_pos, site_size = zip(*sorted_sites)
                site = list(site)
                site_pos = list(site_pos)
                site_size = list(site_size)
                # Add boundaries to the list
                first_pos = site_pos[0]
                first_size = site_size[0]
                last_pos = site_pos[-1]
                last_size = site_size[-1]
                site = [0 - first_pos] + site + [len(record_seq) + 1 + last_pos - last_size]
                site_pos = [first_pos] + list(site_pos) + [last_pos]
                site_size = [first_size] + list(site_size) + [last_size]

                # For each site, compute fragment boundaries and GC content
                for i in range(1, len(site) - 1):
                    count = (i - 1) * 2 + 1
                    # For the negative strand fragment: use the previous site as boundary
                    frag_start = site[i - 1] + site_pos[i - 1]
                    frag_end = site[i] + site_size[i] - site_pos[i] - 1
                    frag_len = frag_end - frag_start
                    # Compute GC content for window up to 200 bp before frag_end
                    gc_seq = record_seq[max(frag_end - 200, 0):frag_end]
                    frag_gc = gc_fraction(gc_seq)
                    f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC:.3f}\n".format(
                        num=count, strand="-", chr=seq_record.id, pos=frag_end, fraglen=frag_len, GC=frag_gc))
                    # For the positive strand fragment: use the current site as boundary
                    frag_start = site[i] + site_pos[i]
                    frag_end = site[i + 1] + site_size[i + 1] - site_pos[i + 1] - 1
                    frag_len = frag_end - frag_start
                    gc_seq = record_seq[frag_start:frag_start + 200]
                    frag_gc = gc_fraction(gc_seq)
                    f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC:.3f}\n".format(
                        num=count + 1, strand="+", chr=seq_record.id, pos=frag_start, fraglen=frag_len, GC=frag_gc))
            else:  # For mnase-based processing
                for i, frag_start in enumerate(range(1, len(record_seq) + 1, binsize)):
                    count = i * 2 + 1
                    frag_end = min(frag_start + binsize - 1, len(record_seq))
                    frag_len = binsize - 1
                    frag_gc = gc_fraction(record_seq[frag_start:frag_end])
                    f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC:.3f}\n".format(
                        num=count, strand="-", chr=seq_record.id, pos=frag_end, fraglen=frag_len, GC=frag_gc))
                    f.write("{num}\t{strand}\t{chr}\t{pos}\t{fraglen}\t{GC:.3f}\n".format(
                        num=count + 1, strand="+", chr=seq_record.id, pos=frag_start, fraglen=frag_len, GC=frag_gc))
    return

def main():
    parser = argparse.ArgumentParser(
        description='Generate genomic features using a FASTA file and a restriction enzyme cut site')
    parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="Input FASTA file")
    parser.add_argument("-s", "--seq", dest="seq", required=True, help="RE cut sequence (or comma-separated for multiple enzymes)")
    parser.add_argument("-o", "--out", dest="outfile", required=True, help="Output file")
    parser.add_argument("-p", "--pos", dest="pos", required=True, help="RE cut position")
    parser.add_argument("-b", "--binsize", dest="binsize", required=False, help="Bin size for MNase-based (default: 5Kb)", default="5Kb")
    args = parser.parse_args()
    bin_size = args.binsize.replace('Kb','000')
    bin_size = bin_size.replace('Mb','000000')
    try:
        bin_size = int(bin_size)
    except ValueError:
        sys.exit("Unknown bin size %s, please double check." % args.binsize)
    find_site(args.fasta, args.seq, args.outfile, args.pos, bin_size)

if __name__ == "__main__":
    main()
