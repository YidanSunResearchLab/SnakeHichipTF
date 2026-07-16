#!/usr/bin/env python3
"""
consensus_interactions.py
=========================
Merge significant interaction calls from four HiChIP callers
(FitHiChIP, HiC-DC+, hichipper, MAPS) into a consensus list.

Each caller uses different coordinate conventions. This script snaps every
anchor to the same fixed-bin grid (default 5 kb), then keeps interactions
supported by >= min_support callers.

For each caller, the best-scoring site (lowest FDR) is reported when the
same binned interaction appears multiple times in one caller's output.

Output columns
--------------
chr start1 end1 start2 end2 support callers
fithichip_counts fithichip_fdr
hicdcplus_counts hicdcplus_fdr
hichipper_counts hichipper_fdr
maps_counts maps_fdr
"""

import argparse
import sys
import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def snap(pos, binsize):
    """Snap a genomic position to the left edge of its bin."""
    return (int(pos) // binsize) * binsize


def make_key(chrom, pos1, pos2, binsize):
    """Canonical (ordered) anchor-pair key on the bin grid."""
    b1 = snap(pos1, binsize)
    b2 = snap(pos2, binsize)
    if b1 > b2:
        b1, b2 = b2, b1
    return (str(chrom), b1, b2)


# --------------------------------------------------------------------------- #
# Per-caller loaders
# Each returns a dict: key -> {"counts": int, "fdr": float}
# --------------------------------------------------------------------------- #
def load_fithichip(path, binsize):
    """
    Columns (0-based): chr1=0 s1=1 chr2=3 s2=4 cc=6 Q-Value_Bias=25
    FDR = Q-Value_Bias (bias-corrected BH-adjusted p).
    When the same binned pair appears more than once, keep the row with
    the lowest FDR (most significant).
    """
    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    records = {}
    for _, r in df.iterrows():
        key = make_key(r.iloc[0], r.iloc[1], r.iloc[4], binsize)
        counts = int(r.iloc[6])
        fdr    = float(r.iloc[25]) if len(r) > 25 else np.nan
        if key not in records or fdr < records[key]["fdr"]:
            records[key] = {"counts": counts, "fdr": fdr}
    return records


def load_hicdcplus(path, binsize):
    """
    Columns: chrI=0 startI=1 chrJ=3 startJ=4 counts=7 qvalue=9
    HiC-DC+ start is 0-based (off by 1 from the bin edge).
    """
    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    records = {}
    for _, r in df.iterrows():
        key = make_key(r.iloc[0], int(r.iloc[1]) + 1,
                                   int(r.iloc[4]) + 1, binsize)
        counts = int(r.iloc[7])
        fdr    = float(r.iloc[9])
        if key not in records or fdr < records[key]["fdr"]:
            records[key] = {"counts": counts, "fdr": fdr}
    return records


def load_hichipper(path, binsize):
    """
    hichipper significant interactions BEDPE.
    Columns: chr_1=0 start_1=1 end_1=2 chr_2=3 start_2=4 end_2=5
             sample_name=6  mango.FDR=7
    Note: the significant.bedpe does not carry raw contact counts;
    the sample_name string is in col 6. FDR = mango.FDR (col 7).
    Anchors are variable-width peaks; midpoint is snapped to the bin grid.
    """
    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    records = {}
    for _, r in df.iterrows():
        mid1 = (int(r.iloc[1]) + int(r.iloc[2])) // 2
        mid2 = (int(r.iloc[4]) + int(r.iloc[5])) // 2
        key  = make_key(r.iloc[0], mid1, mid2, binsize)
        # col 6 is a string (sample name), col 7 is mango FDR
        fdr = float(r.iloc[7]) if len(r) > 7 else np.nan
        if key not in records or (not np.isnan(fdr) and
                                  (np.isnan(records[key]["fdr"]) or
                                   fdr < records[key]["fdr"])):
            records[key] = {"counts": "none", "fdr": fdr}
    return records


def load_maps(path, binsize):
    """
    Columns: chr1=0 start1=1 chr2=3 start2=4 count=6 fdr=8
    Anchors are already on the bin grid.
    """
    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    records = {}
    for _, r in df.iterrows():
        key = make_key(r.iloc[0], r.iloc[1], r.iloc[4], binsize)
        counts = int(r.iloc[6])
        fdr    = float(r.iloc[8])
        if key not in records or fdr < records[key]["fdr"]:
            records[key] = {"counts": counts, "fdr": fdr}
    return records


# --------------------------------------------------------------------------- #
LOADERS = {
    "fithichip":  load_fithichip,
    "hicdcplus":  load_hicdcplus,
    "hichipper":  load_hichipper,
    "maps":       load_maps,
}

CALLER_DISPLAY = {
    "fithichip": "FitHiChIP",
    "hicdcplus": "HiC-DC+",
    "hichipper": "hichipper",
    "maps":      "MAPS",
}


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fithichip",   default=None)
    ap.add_argument("--hicdcplus",  default=None)
    ap.add_argument("--hichipper",  default=None)
    ap.add_argument("--maps",       default=None)
    ap.add_argument("--binsize",    type=int, default=5000)
    ap.add_argument("--min-support",type=int, default=2,
                    help="minimum callers required [2]")
    ap.add_argument("--out",        required=True)
    args = ap.parse_args()

    # ---- load each caller ----
    caller_data = {}   # name -> {key: {counts, fdr}}
    for name, loader in LOADERS.items():
        path = getattr(args, name)
        if not path:
            continue
        try:
            records = loader(path, args.binsize)
            caller_data[name] = records
            print(f"[consensus] {CALLER_DISPLAY[name]}: "
                  f"{len(records)} interactions", file=sys.stderr)
        except Exception as e:
            print(f"[consensus] {CALLER_DISPLAY[name]}: FAILED ({e})",
                  file=sys.stderr)

    if not caller_data:
        sys.exit("[consensus] ERROR: no caller files loaded")

    # ---- aggregate across callers ----
    all_keys = set().union(*[set(d) for d in caller_data.values()])
    print(f"[consensus] union: {len(all_keys)} interactions", file=sys.stderr)

    rows = []
    for key in all_keys:
        chrom, start1, start2 = key
        supporting = [c for c, d in caller_data.items() if key in d]
        n = len(supporting)

        row = {
            "chr":     chrom,
            "start1":  start1,
            "end1":    start1 + args.binsize,
            "start2":  start2,
            "end2":    start2 + args.binsize,
            "support": n,
            "callers": ",".join(CALLER_DISPLAY[c] for c in sorted(supporting)),
        }
        # per-caller stats
        for name in LOADERS:
            if name in caller_data and key in caller_data[name]:
                info = caller_data[name][key]
                row[f"{name}_counts"] = info["counts"]
                row[f"{name}_fdr"]    = info["fdr"]
            else:
                row[f"{name}_counts"] = "none"
                row[f"{name}_fdr"]    = "none"
        rows.append(row)

    col_order = ["chr", "start1", "end1", "start2", "end2",
                 "support", "callers",
                 "fithichip_counts", "fithichip_fdr",
                 "hicdcplus_counts",  "hicdcplus_fdr",
                 "hichipper_counts", "hichipper_fdr",
                 "maps_counts",      "maps_fdr"]

    df = (pd.DataFrame(rows)[col_order]
            .sort_values(["chr", "start1", "start2"]))

    # full table
    all_out = args.out.replace(".bedpe", ".all_callers.bedpe")
    df.to_csv(all_out, sep="\t", index=False)

    # consensus table
    consensus = df[df["support"] >= args.min_support]
    consensus.to_csv(args.out, sep="\t", index=False)

    print(f"[consensus] consensus (>={args.min_support}): "
          f"{len(consensus)} -> {args.out}", file=sys.stderr)
    print("\n[consensus] support distribution:", file=sys.stderr)
    print(df["support"].value_counts().sort_index().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()