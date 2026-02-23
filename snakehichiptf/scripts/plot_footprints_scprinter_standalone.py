import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyBigWig
from pathlib import Path

# ================= PATHS =================
main_dir = Path("/home/syidan/syidan/Data/Processed/scPrinter/brain_difPart/Caudate")
work_dir = main_dir / "seq2print"
plots_dir = work_dir / "plots"
regions_bed = main_dir / "RHOB_regions.bed"
# ========================================

os.makedirs(plots_dir, exist_ok=True)

# ---------------- Samples ----------------
all_samples = sorted([
    d.name.replace("tfbs_", "")
    for d in work_dir.iterdir()
    if d.is_dir() and d.name.startswith("tfbs_")
])

if len(all_samples) == 0:
    raise RuntimeError("No tfbs_* directories found")

# Use only the first sample
sample = all_samples[0]
samples = [sample]
print("Using first sample only:", sample)

# ---------------- Regions ----------------
regions = pd.read_csv(
    regions_bed,
    sep="\t",
    header=None,
    names=["chrom", "start", "end"]
)
regions["name"] = [f"RHOB_{i+1}" for i in range(len(regions))]

# ---------------- Heatmaps ----------------
for _, r in regions.iterrows():
    chrom, start, end, name = r.chrom, r.start, r.end, r.name
    print(f"Processing {name}: {chrom}:{start}-{end}")

    tfbs_matrix = []

    for sample in samples:
        bw_path = work_dir / f"tfbs_{sample}" / f"{sample}_TFBS.bigwig"

        with pyBigWig.open(str(bw_path)) as bw:
            # Handle chr / no-chr automatically
            chrom_use = chrom
            if chrom_use not in bw.chroms():
                chrom_use = chrom.replace("chr", "") if chrom.startswith("chr") else f"chr{chrom}"

            signal = np.array(
                bw.values(chrom_use, start, end, numpy=True)
            ).squeeze()  # ensure 1D

        tfbs_matrix.append(signal)

    # Build (samples × positions) matrix
    tfbs = pd.DataFrame(
        np.vstack(tfbs_matrix),
        index=samples,
        columns=np.arange(start, end)
    )

    # ---------------- COLOR SCALING ----------------
    vmax = np.nanpercentile(tfbs.values, 99)

    # ---------------- Plot ----------------
    plt.figure(figsize=(4, 2))
    sns.heatmap(
        tfbs,
        cmap="RdBu_r",
        vmin=0,
        vmax=vmax,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={"label": "Predicted TFBS signal"}
    )

    plt.axis('off')       # remove axes
    plt.tight_layout()

    # Save as PDF
    out = plots_dir / f"TFBS_{name}_{sample}.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()

    print(f"  Saved {out.name}\n")

print("✅ Done — first-sample TFBS heatmaps generated as PDF.")
