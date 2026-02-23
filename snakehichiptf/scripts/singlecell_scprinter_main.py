import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"  # default: use one GPU; change if needed

import torch
import scprinter as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import numpy as np
import pickle
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from scanpy.plotting.palettes import zeileis_28
import anndata
import scanpy as sc
import json
import csv
import re
from pathlib import Path
import argparse

import snapatac2 as snap
import pyranges as pr
from tqdm.auto import tqdm

import pyBigWig as pw
from sklearn.neighbors import NearestNeighbors

import glob

print("[DEBUG] >>> USING NEW SINGLECELL SCRIPT V2024-12-05B <<<", flush=True)

# -------------------------------------------------------------------
# Patch SnapATAC2 TSSe to avoid Rust panic on tiny / test data
# -------------------------------------------------------------------
def _dummy_tsse(adata, gene_anno):
    """
    Dummy TSSe function: always returns zeros.

    This bypasses the SnapATAC2 Rust implementation that may panic
    with 'Unknown frame descriptor' on very small datasets.
    """
    return np.zeros(adata.n_obs, dtype=float)


snap.metrics.tsse = _dummy_tsse


# -------------------------------------------------------------------
# Simple file checkpoint helper
# -------------------------------------------------------------------
def exists_nonempty(path: str) -> bool:
    """
    Return True if file exists and has non-zero size.
    Used as a lightweight checkpoint to skip finished steps.
    """
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0


# -------------------------------------------------------------------
# Fetch bigwig signal for TFBS sites
# -------------------------------------------------------------------
def fetch_bw(args):
    """
    Given (TFBS dataframe, bigwig path, genome), return the mean signal
    in each TFBS interval.
    """
    TFBS, bw, genome = args
    chroms, starts, ends = np.array(TFBS["chrom"]), np.array(TFBS["start"]), np.array(TFBS["end"])
    res_all = {}

    with pw.open(bw, "r") as f:
        for chrom in tqdm(genome.chrom_sizes):
            if chrom == "chrY":
                # skip chrY to avoid issues when missing
                continue
            res_all[chrom] = f.values(chrom, 0, genome.chrom_sizes[chrom], numpy=True)

    vs = []
    for chrom, left, right in zip(tqdm(chroms, mininterval=1), starts, ends):
        vs.append(np.nanmean(res_all[chrom][left:right]))
    return vs


# -------------------------------------------------------------------
# Select the first N peak regions that overlap footprint intervals
# -------------------------------------------------------------------
def select_regions_with_footprints(peaks_bed, footprints_bed, out_bed, n_regions=10):
    """
    Intersect peaks with footprints and write the first N intervals
    to out_bed.
    """
    peaks = pd.read_csv(
        peaks_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )
    footprints = pd.read_csv(
        footprints_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )

    peaks_pr = pr.PyRanges(peaks)
    fp_pr = pr.PyRanges(footprints)

    overlapped = peaks_pr.join(fp_pr)
    if len(overlapped) == 0:
        raise ValueError(
            f"No overlap between {peaks_bed} and {footprints_bed}."
        )

    df = overlapped.df[["Chromosome", "Start", "End"]].drop_duplicates()
    df = df.sort_values(["Chromosome", "Start", "End"])
    selected = df.head(n_regions)

    selected.to_csv(out_bed, sep="\t", header=False, index=False)


# -------------------------------------------------------------------
# Construct pseudobulk groups
# -------------------------------------------------------------------

def build_pseudobulk_groups(chromvar, adata, depth_threshold=5e6, n_centers_per_type=10):
    print(f"[INFO] Building pseudobulk groups "
          f"(depth_threshold={depth_threshold}, centers_per_type={n_centers_per_type}) ...")

    # Fragment depth per cell, based on the cell-by-peak matrix
    cell_depth = np.array(adata.X.sum(axis=1)).squeeze()
    cell_barcodes = np.array(chromvar.obs.index)

    # Fit KNN in chromVAR PCA space
    if "X_pca" not in chromvar.obsm:
        raise ValueError("chromvar.obsm['X_pca'] is missing. Run PCA before building pseudobulks.")
    nbrs = NearestNeighbors(
        n_neighbors=chromvar.shape[0],
        algorithm="ball_tree"
    ).fit(chromvar.obsm["X_pca"])

    pbulk_centers = []
    barcode_groups = []

    celltypes = np.unique(chromvar.obs["celltype"])
    for celltype in celltypes:
        celltype_center_count = 0
        cell_inds = np.where(chromvar.obs["celltype"] == celltype)[0]
        if len(cell_inds) == 0:
            continue

        print(f"[INFO]   celltype={celltype}, n_cells={len(cell_inds)}")

        # Keep sampling centers until we have enough pseudobulks for this cell type
        while celltype_center_count < n_centers_per_type:
            # Avoid infinite loops if this cell type is very rare
            if len(cell_inds) < 5:
                print(f"[WARN]   celltype={celltype} has too few cells; "
                      f"stop at {celltype_center_count} pseudobulks.")
                break

            # 1) Randomly choose a center cell of this cell type
            new_center = np.random.choice(cell_inds, 1)[0]

            # 2) Order all cells by distance to the center in PCA space
            distances, indices = nbrs.kneighbors(
                [chromvar.obsm["X_pca"][new_center, :]]
            )
            nbr_inds = indices[0, :]

            # 3) Accumulate fragment depth until we reach the threshold
            cumulative_depth = np.cumsum(cell_depth[nbr_inds])

            # If total depth is still below threshold, use all neighbors
            if cumulative_depth[-1] <= depth_threshold:
                n_members = len(nbr_inds)
            else:
                # First index where cumulative depth exceeds the threshold
                n_members = int(np.where(cumulative_depth > depth_threshold)[0][0])

            pbulk_members = nbr_inds[:n_members]

            # 4) Compute purity of this pseudobulk
            celltype_labels = chromvar.obs["celltype"].values[pbulk_members]
            vals, counts = np.unique(celltype_labels, return_counts=True)
            majority_type = vals[np.argmax(counts)]
            purity = np.max(counts) / len(celltype_labels)

            # Discard impure pseudobulks
            if purity < 0.99:
                continue

            # 5) Accept this pseudobulk
            pbulk_centers.append(new_center)
            celltype_center_count += 1

            center_ind = len(pbulk_centers)
            chromvar.obs[f"pbulk_{center_ind}_member"] = np.zeros(chromvar.shape[0])
            chromvar.obs[f"pbulk_{center_ind}_member"][pbulk_members] = 1

            # Record barcode-to-group mapping
            new_bc_group = pd.DataFrame({
                "barcode": cell_barcodes[pbulk_members],
                "group": f"{celltype}_pbulk_{celltype_center_count}",
            })
            barcode_groups.append(new_bc_group)

        print(f"[INFO]   built {celltype_center_count} pseudobulks for {celltype}")

    if len(barcode_groups) == 0:
        raise RuntimeError(
            "[ERROR] No pseudobulk groups were constructed; "
            "check cell types and depth_threshold."
        )

    barcode_groups = pd.concat(barcode_groups, axis=0)
    return barcode_groups, np.array(pbulk_centers)


## -------------------------------------------------------------------
# Train base + LoRA seq2PRINT models locally (no wandb)
## -------------------------------------------------------------------
def train_base_and_lora_locally(
    printer,
    genome,
    chromvar,
    cell_grouping,        # pseudobulk grouping from df2cell_grouping
    group_names,
    work_dir,
    cleaned_seq2print_bed,
    gpu_id=0,
):
    """
    Train a base seq2PRINT model (single bulk group) and then a LoRA model
    on top of it using pseudobulk groupings, all locally without wandb.

    Now added:
    - If a base model already exists in work_dir/model, skip retraining base;
    - If a LoRA model already exists, skip retraining LoRA.
    """

    configs_dir = os.path.join(work_dir, "configs")
    model_dir = os.path.join(work_dir, "model")
    temp_dir = os.path.join(work_dir, "temp")
    os.makedirs(configs_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    # ---------------------------------------------------
    # 1) Base model (all cells as a single bulk group)
    # ---------------------------------------------------
    base_cfg_path = os.path.join(configs_dir, "PBMC_Bulk_local_fold0.json")
    print("[INFO] Creating base seq2PRINT config:", base_cfg_path)

    base_cell_grouping = printer.obs_names   # all barcodes
    base_group_name = "brain_Bulk"          # just a single group name (string)

    base_config = scp.tl.seq_model_config(
        printer=printer,
        region_path=cleaned_seq2print_bed,
        cell_grouping=base_cell_grouping,   # single group
        group_names=base_group_name,        # string, not a list
        genome=genome,
        fold=0,
        overwrite_bigwig=False,
        model_name="PBMC_Bulk_local",
        additional_config={
            "notes": "local_training_no_wandb",
            "tags": ["PBMC", "Bulk", "local"],
        },
        path_swap=(work_dir, ""),
        config_save_path=base_cfg_path,
    )

    # If a base model already exists, skip retraining
    base_candidates = [
        os.path.join(model_dir, f)
        for f in os.listdir(model_dir)
        if f.endswith(".pt") and "LoRA" not in f
    ]
    if base_candidates:
        print("[INFO] Found existing base model checkpoint(s); skip base training.")
    else:
        print("[INFO] Training base seq2PRINT model (local, no wandb) ...")
        scp.tl.launch_seq2print(
            model_config_path=base_cfg_path,
            temp_dir=temp_dir,
            model_dir=model_dir,
            data_dir=work_dir,
            gpus=gpu_id,
            wandb_project=None,  # disable wandb
            verbose=True,
            launch=True,
        )
        base_candidates = [
            os.path.join(model_dir, f)
            for f in os.listdir(model_dir)
            if f.endswith(".pt") and "LoRA" not in f
        ]

    if not base_candidates:
        raise RuntimeError("[ERROR] No .pt base model found in model_dir after training.")
    base_candidates.sort(key=os.path.getmtime)
    base_model_path = base_candidates[-1]
    print("[INFO] Using base model checkpoint:", base_model_path)

    # ---------------------------------------------------
    # 2) LoRA config: Use pseudobulk + PCA embeddings
    # ---------------------------------------------------
    embeddings = pd.DataFrame(chromvar.obsm["X_pca"], index=chromvar.obs_names)
    lora_cfg_path = os.path.join(configs_dir, "PBMC_LoRA_local_fold0.json")

    print("[INFO] Creating LoRA config:", lora_cfg_path)
    lora_config = scp.tl.seq_lora_model_config(
        printer=printer,
        region_path=cleaned_seq2print_bed,
        cell_grouping=cell_grouping,        # pseudobulk grouping
        group_names=group_names,            # multiple pseudobulk group names
        embeddings=embeddings,
        genome=genome,
        pretrain_model=base_model_path,
        overwrite_barcode=False,
        model_name="PBMC_LoRA_local",
        fold=0,
        model_config=base_cfg_path,
        additional_lora_config={
            "lr": 3e-5,
            "notes": "local_lora_no_wandb",
            "tags": ["PBMC", "LoRA", "local"],
        },
        path_swap=(work_dir, ""),
        config_save_path=lora_cfg_path,
    )

    lora_candidates = [
        os.path.join(model_dir, f)
        for f in os.listdir(model_dir)
        if f.endswith(".pt") and "LoRA" in f
    ]
    if lora_candidates:
        print("[INFO] Found existing LoRA checkpoint(s); skip LoRA training.")
    else:
        print("[INFO] Training LoRA model (local, no wandb) ...")
        scp.tl.launch_seq2print(
            model_config_path=lora_cfg_path,
            temp_dir=temp_dir,
            model_dir=model_dir,
            data_dir=work_dir,
            gpus=gpu_id,
            wandb_project=None,
            verbose=True,
            launch=True,
        )
        lora_candidates = [
            os.path.join(model_dir, f)
            for f in os.listdir(model_dir)
            if f.endswith(".pt") and "LoRA" in f
        ]

    if not lora_candidates:
        print("[WARN] No LoRA-specific checkpoint found; fall back to newest base model.")
        lora_candidates = base_candidates

    lora_candidates.sort(key=os.path.getmtime)
    lora_model_path = lora_candidates[-1]
    print("[INFO] Using LoRA model checkpoint:", lora_model_path)

    return lora_model_path, lora_cfg_path


# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Single-cell scPrinter pipeline (PBMC tutorial style, end-to-end)."
    )
    p.add_argument(
        "--main_dir",
        required=True,
        help="Project main directory; will contain seq2print/ and fragments/."
    )
    p.add_argument(
        "--fragment_file",
        required=True,
        help="Single-cell fragments.tsv.gz (columns: chrom, start, end, barcode)."
    )
    p.add_argument(
        "--device",
        default="cuda",
        choices=["cuda", "cpu"],
        help="Device for chromVAR deviations."
    )
    p.add_argument(
        "--min_fragments",
        type=int,
        default=1000,
        help="Minimum number of fragments per cell in import_fragments."
    )
    p.add_argument(
        "--min_tsse",
        type=float,
        default=7.0,
        help="Minimum TSSe in import_fragments."
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="If set, skip steps whose outputs already exist."
    )
    # --- LoRA control (local model, no wandb) ---
    p.add_argument(
        "--no_lora",
        action="store_true",
        help="If set, skip LoRA + TFBS pseudobulk part even if model/config are provided."
    )
    p.add_argument(
        "--lora_model_path",
        default=None,
        help="Path to local LoRA model (.pt). If None, LoRA is skipped unless overridden."
    )
    p.add_argument(
        "--lora_config_path",
        default=None,
        help="Path to local LoRA config JSON. If None, LoRA is skipped unless overridden."
    )
    p.add_argument(
        "--regions_bed",
        default=None,
        help="Optional BED with regions of interest for LoRA TFBS visualization."
    )
    p.add_argument(
        "--footprints_bed",
        default=None,
        help="Optional BED with footprints; if provided, select regions that overlap footprints."
    )
    p.add_argument(
        "--n_regions",
        type=int,
        default=10,
        help="Number of regions to select automatically when regions_bed is not given."
    )
    p.add_argument(
        "--marker_dict",
        type=str,
        default="{}",
        help="JSON-style dict mapping TF motifs to cell types"
    )


    return p.parse_args()


# -------------------------------------------------------------------
# Main pipeline
# -------------------------------------------------------------------
def main():
    args = parse_args()

    main_dir = args.main_dir.rstrip("/")
    fragment_file = os.path.abspath(args.fragment_file)
    device = args.device

    work_dir = f"{main_dir}/seq2print"
    frag_dir = f"{main_dir}/fragments"
    plots_dir = f"{work_dir}/plots"
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(frag_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    configs_dir = os.path.join(work_dir, "configs")
    model_dir = os.path.join(work_dir, "model")
    os.makedirs(configs_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)

    # Symlink the single-cell fragments into fragments/
    frag_basename = os.path.basename(fragment_file)
    frag_link = os.path.join(frag_dir, frag_basename)
    if not os.path.exists(frag_link):
        os.symlink(fragment_file, frag_link)

    genome = scp.genome.hg38

    # ------------------------------------------------------------
    # 1. Import fragments into scPrinter object
    # ------------------------------------------------------------
    printer_path = os.path.join(work_dir, "PBMC_scATAC_scprinter.h5ad")
    if args.resume and exists_nonempty(printer_path):
        print(f"[INFO] {printer_path} exists, loading printer (resume).")
        printer = scp.load_printer(printer_path, genome)
    else:
        print("[INFO] Importing single-cell fragments into scPrinter object ...")
        printer = scp.pp.import_fragments(
            path_to_frags=fragment_file,
            barcodes=None,
            savename=printer_path,
            genome=genome,
            min_num_fragments=args.min_fragments,
            min_tsse=args.min_tsse,
            sorted_by_barcode=False,
            low_memory=False,
        )

    # ------------------------------------------------------------
    # 2. Peak calling (seq2PRINT + chromvar presets)
    # ------------------------------------------------------------
    cleaned_seq2print_bed = os.path.join(work_dir, "seq2print_cleaned_narrowPeak.bed")
    chromvar_regions_bed = os.path.join(work_dir, "regions.bed")

    # 2.1 preset='seq2PRINT'
    if args.resume and exists_nonempty(cleaned_seq2print_bed):
        print(f"[INFO] {cleaned_seq2print_bed} exists, skip seq2PRINT peak calling.")
    else:
        print("[INFO] Calling peaks with preset='seq2PRINT' ...")
        scp.pp.call_peaks(
            printer=printer,
            frag_file=fragment_file,
            cell_grouping=[None],
            group_names=["all"],
            preset="seq2PRINT",
            n_jobs=1,
        )
        cleaned_peaks = pd.DataFrame(printer.uns["peak_calling"]["all_cleaned"][:])
        cleaned_peaks.to_csv(
            cleaned_seq2print_bed,
            sep="\t",
            header=False,
            index=False,
        )

    # 2.2 preset='chromvar'
    if args.resume and exists_nonempty(chromvar_regions_bed):
        print(f"[INFO] {chromvar_regions_bed} exists, skip chromvar peak calling.")
    else:
        print("[INFO] Calling peaks with preset='chromvar' ...")
        scp.pp.call_peaks(
            printer=printer,
            frag_file=fragment_file,
            cell_grouping=[None],
            group_names=["chromvar_all"],
            preset="chromvar",
            overwrite=False,
        )
        chromvar_peaks = pd.DataFrame(
            printer.uns["peak_calling"]["chromvar_all_cleaned"][:]
        )
        chromvar_peaks.to_csv(
            chromvar_regions_bed,
            sep="\t",
            header=False,
            index=False,
        )

    # ------------------------------------------------------------
    # 3. Cell-by-peak matrix and basic coverage filtering
    # ------------------------------------------------------------
    cell_peak_path = os.path.join(work_dir, "cell_peak.h5ad")
    if args.resume and exists_nonempty(cell_peak_path):
        print(f"[INFO] {cell_peak_path} exists, skip make_peak_matrix.")
    else:
        print("[INFO] Building cell-by-peak matrix from regions.bed ...")
        adata = scp.pp.make_peak_matrix(
            printer,
            regions=chromvar_regions_bed,
            region_width=300,
            cell_grouping=None,
            group_names=None,
            sparse=True,
        )
        adata.write(cell_peak_path)

    print("[INFO] Filtering peaks with zero coverage ...")
    adata = anndata.read_h5ad(cell_peak_path)
    coverage = np.array(adata.X.sum(axis=0)).squeeze()
    adata = adata[:, coverage > 0]
    adata.write(cell_peak_path)

    # ------------------------------------------------------------
    # Save filtered regions to regions_filt.bed (used by Snakemake)
    # ------------------------------------------------------------
    regions = pd.read_csv(chromvar_regions_bed, sep="\t", header=None)
    peak_depth = np.array(adata.X.sum(axis=0)).squeeze()

    # You can tune the threshold; 200 follows the PBMC tutorial pattern
    keep_idx = np.where(peak_depth > 200)[0]
    regions_filt = regions.iloc[keep_idx, :]

    regions_filt_bed = os.path.join(work_dir, "regions_filt.bed")
    regions_filt.to_csv(
        regions_filt_bed,
        sep="\t",
        header=False,
        index=False,
    )

    # ------------------------------------------------------------
    # 4. chromVAR motif deviations + clustering + cell-type labels
    # ------------------------------------------------------------
    chromvar_h5ad = os.path.join(work_dir, "chromvar_cisbp.h5ad")
    if args.resume and exists_nonempty(chromvar_h5ad):
        print(f"[INFO] {chromvar_h5ad} exists, loading chromVAR result.")
        chromvar = anndata.read_h5ad(chromvar_h5ad)
    else:
        print("[INFO] Running chromVAR background sampling ...")
        scp.chromvar.sample_bg_peaks(
            adata,
            genome=genome,
            method="chromvar",
            niterations=250,
        )

        print("[INFO] Scanning motifs (FigR Human Motifs) ...")
        motif = scp.motifs.FigR_Human_Motifs(
            genome,
            bg=list(adata.uns["bg_freq"]),
            n_jobs=32,
            pvalue=5e-5,
            mode="motifmatchr",
        )
        motif.prep_scanner(None, pvalue=5e-5)
        motif.chromvar_scan(adata)

        print(f"[INFO] Computing chromVAR deviations on device={device} ...")
        chromvar = scp.chromvar.compute_deviations(
            adata,
            chunk_size=50000,
            device=device,
        )
        chromvar.write(chromvar_h5ad)

    print(f"[INFO] chromVAR deviations matrix shape: {chromvar.X.shape}")

    # PCA / neighbors / UMAP / Leiden clustering
    print("[INFO] Running PCA / neighbors / UMAP / Leiden ...")
    sc.tl.pca(chromvar)
    sc.pp.neighbors(chromvar)
    sc.tl.umap(chromvar)
    sc.tl.leiden(
        chromvar,
        flavor="igraph",
        n_iterations=10,
        resolution=0.2,
        random_state=2,
    )

    # UMAP colored by selected TF motifs
    for col in ["CEBPA", "TCF7", "PAX5"]:
        if col in chromvar.var_names:
            sc.pl.umap(
                chromvar,
                color=col,
                cmap="RdBu_r",
                vmin=-3,
                vmax=3,
                size=2,
                show=False,
             )
            plt.savefig(os.path.join(plots_dir, f"umap_{col}.png"), dpi=200)
            plt.close()

    sc.pl.umap(chromvar, color="leiden", size=2, show=False)
    plt.savefig(os.path.join(plots_dir, "umap_leiden.png"), dpi=200)
    plt.close()

    # UMAP
    raw = args.marker_dict
    if raw is None:
    	marker_dict = {}
    else:
    	marker_dict = json.loads(raw)
    print("[INFO] marker_dict loaded:", marker_dict)

    if len(marker_dict) == 0:
        print("[WARN] marker_dict is empty. Falling back to single label assignment.")
        chromvar.obs["celltype"] = "unknown"
        chromvar.obs["celltype"] = chromvar.obs["celltype"].astype("category")

    else:
        print("[INFO] Assigning cell types from marker_dict...")

        marker_TFs = list(marker_dict.keys())
        clusters = np.unique(chromvar.obs["leiden"])

        cluster_score = []

        for cluster in clusters:
            cell_idx = chromvar.obs["leiden"] == cluster
            cluster_score.append(np.mean(chromvar[cell_idx, marker_TFs].X, axis=0))

        cluster_score = np.array(cluster_score)

        map_dict = {}
        for i, cluster in enumerate(clusters):
            best_tf = marker_TFs[np.argmax(cluster_score[i, :])]
            map_dict[cluster] = marker_dict[best_tf]

        print("[INFO] Cluster → celltype mapping:", map_dict)
        chromvar.obs["celltype"] = [
            map_dict[cl] for cl in chromvar.obs["leiden"]
        ]
        chromvar.obs["celltype"] = chromvar.obs["celltype"].astype("category")

    sc.pl.umap(chromvar, color="celltype", size=2, show=False)
    plt.savefig(os.path.join(plots_dir, "umap_celltype.png"), dpi=200)
    plt.close()
    chromvar.write(chromvar_h5ad)
    # ------------------------------------------------------------
    # 5.2 Build pseudobulk groups (official-like, but without wandb)
    # ------------------------------------------------------------
    print("[INFO] Constructing pseudobulk groups (official-style) ...")

    adata = anndata.read_h5ad(cell_peak_path)

    barcode_groups_path = os.path.join(work_dir, "barcode_groups.csv")
    pbulk_centers_path = os.path.join(work_dir, "pbulk_centers.npy")

    if args.resume and exists_nonempty(barcode_groups_path) and exists_nonempty(pbulk_centers_path):
        print(f"[INFO] Found existing pseudobulk groups, loading from disk:")
        print(f"       {barcode_groups_path}")
        print(f"       {pbulk_centers_path}")
        barcode_groups = pd.read_csv(barcode_groups_path)
        pbulk_centers = np.load(pbulk_centers_path)
    else:
        barcode_groups, pbulk_centers = build_pseudobulk_groups(
            chromvar=chromvar,
            adata=adata,
            depth_threshold=5e6,       # change to ~5e3 for toy tests
            n_centers_per_type=10,
        )
        # Save checkpoint
        barcode_groups.to_csv(barcode_groups_path, index=False)
        np.save(pbulk_centers_path, pbulk_centers)

    # Visualize pseudobulk centers on UMAP (black dots = pseudobulk centers)
    chromvar.obs["pbulk_center"] = np.zeros(chromvar.shape[0])
    chromvar.obs["pbulk_center"][pbulk_centers] = 1
    import matplotlib.colors as mcolors
    cmap = mcolors.ListedColormap(["#D3D3D3", "black"])
    sc.pl.umap(chromvar, color="pbulk_center", size=5, cmap=cmap, show=False)
    plt.savefig(os.path.join(plots_dir, "umap_pbulk_centers.png"), dpi=200)
    plt.close()

    # Convert barcode_groups to the format required by scPrinter
    cell_grouping, group_names = scp.utils.df2cell_grouping(printer, barcode_groups)

    # ------------------------------------------------------------
    # 5.3 Decide how to obtain a LoRA model:
    #     - If user provides --lora_model_path / --lora_config_path and they exist,
    #       use them directly.
    #     - Otherwise, train a base seq2PRINT model and a LoRA model locally.
    # ------------------------------------------------------------
    lora_model_path = args.lora_model_path
    lora_config_path = args.lora_config_path

    lora_model_path = args.lora_model_path
    lora_config_path = args.lora_config_path

    # If user did not explicitly provide, and there is already a LoRA model and default config in work_dir, reuse them directly
    if not lora_model_path and not lora_config_path:
        default_lora_cfg = os.path.join(configs_dir, "PBMC_LoRA_local_fold0.json")
        existing_lora_models = []
        if os.path.exists(model_dir):
            existing_lora_models = [
                os.path.join(model_dir, f)
                for f in os.listdir(model_dir)
                if f.endswith(".pt") and "LoRA" in f
            ]
        if existing_lora_models and exists_nonempty(default_lora_cfg):
            existing_lora_models.sort(key=os.path.getmtime)
            lora_model_path = existing_lora_models[-1]
            lora_config_path = default_lora_cfg
            print("[INFO] Found existing LoRA model/config in work_dir; will reuse:")
            print("       model :", lora_model_path)
            print("       config:", lora_config_path)

    if args.no_lora:
        print("[INFO] --no_lora set; finish after preprocessing and pseudobulk.")
        print(f"[INFO] Outputs under: {work_dir}")
        return

    # Case 1: user-provided LoRA model/config and both paths exist
    if lora_model_path and lora_config_path:
        lora_model_path = os.path.abspath(lora_model_path)
        lora_config_path = os.path.abspath(lora_config_path)

        if os.path.exists(lora_model_path) and os.path.exists(lora_config_path):
            print("[INFO] Using user-provided LoRA model/config:")
            print("       model :", lora_model_path)
            print("       config:", lora_config_path)
        else:
            print("[WARN] Provided LoRA model/config paths do not exist; "
                  "fall back to local training.")
            lora_model_path = None
            lora_config_path = None

    # Case 2: no valid LoRA provided -> train base + LoRA locally
    if not lora_model_path or not lora_config_path:
        print("[INFO] No valid LoRA model/config detected.")
        print("[INFO] Training base seq2PRINT and LoRA locally (no wandb).")

        lora_model_path, lora_config_path = train_base_and_lora_locally(
            printer=printer,
            genome=genome,
            chromvar=chromvar,
            cell_grouping=cell_grouping,
            group_names=group_names,
            work_dir=work_dir,
            cleaned_seq2print_bed=cleaned_seq2print_bed,
            gpu_id=0,   # adjust GPU index if needed
        )

    # Now load the LoRA config
    lora_model_path = os.path.abspath(lora_model_path)
    lora_config_path = os.path.abspath(lora_config_path)

    import json as _json
    with open(lora_config_path) as f:
        lora_cfg = _json.load(f)

    # Ensure the config has group_names consistent with our pseudobulks
    if "group_names" not in lora_cfg or not lora_cfg["group_names"]:
        lora_cfg["group_names"] = list(group_names)

    print("[INFO] Using LoRA model:", lora_model_path)
    print("[INFO] Using LoRA config:", lora_config_path)


    # 5.5 Select regions for TFBS visualization 

    regions_bed_arg = args.regions_bed
    footprints_bed = args.footprints_bed
    n_regions = args.n_regions

    regions_test_bed = os.path.join(work_dir, "regions_test.bed")

    # Base selection on regions_filt.bed (already coverage-filtered)
    cleaned_peaks_bed = cleaned_seq2print_bed
    # or  use regions_filt_bed:
    peaks_source_bed = regions_filt_bed  # or cleaned_peaks_bed

    if regions_bed_arg and regions_bed_arg.strip():
        # 1) User-provided regions_bed
        print(f"[INFO] Using user-specified regions_bed: {regions_bed_arg}")
        regions_df = pd.read_csv(regions_bed_arg, sep="\t", header=None).iloc[:, :3]
        regions_df.columns = ["chrom", "start", "end"]
        regions_df.to_csv(regions_test_bed, sep="\t", header=False, index=False)

    else:
        # 2) Auto-select regions
        if footprints_bed and footprints_bed.strip():
            print(f"[INFO] Using footprints_bed={footprints_bed} to select regions.")
            select_regions_with_footprints(
                peaks_bed=peaks_source_bed,
                footprints_bed=footprints_bed,
                out_bed=regions_test_bed,
                n_regions=n_regions,
            )
            regions_df = pd.read_csv(regions_test_bed, sep="\t", header=None).iloc[:, :3]
        else:
            print("[INFO] No regions_bed or footprints_bed provided; using first N peaks.")
            peaks = pd.read_csv(peaks_source_bed, sep="\t", header=None)
            regions_df = peaks.iloc[:n_regions, :3]
            regions_df.to_csv(
                regions_test_bed,
                sep="\t",
                header=False,
                index=False,
            )

        regions_df.columns = ["chrom", "start", "end"]

    # Build regions_dict in the same style as before, but generic
    regions_dict = {
        f"{row.chrom}:{row.start}-{row.end}": f"region_{i+1}"
        for i, row in regions_df.iterrows()
    }

    regions_df = []
    for region in regions_dict:
        regions_df.append(re.split("[:-]", region))
    regions_df = pd.DataFrame(regions_df)
    regions_test_bed = os.path.join(work_dir, "regions_test.bed")
    regions_df.to_csv(regions_test_bed, sep="\t", header=False, index=False)

    # 5.6 Run LoRA TFBS on the small example regions and plot heatmaps
    print("[INFO] Running LoRA TFBS on regions_test.bed for visualization ...")
    adata_tfbs = scp.tl.seq_tfbs_seq2print(
        seq_attr_count=None,
        seq_attr_footprint=None,
        genome=printer.genome,
        region_path=regions_test_bed,
        gpus=[0],
        model_type="lora",
        model_path=lora_model_path,
        lora_config=lora_cfg,          # patched with group_names above
        group_names=list(group_names),
        verbose=False,
        launch=True,                  # <<< run in-process, no CLI
        return_adata=True,
        overwrite_seqattr=True,
        save_key="PBMC_LoRA_test",
        save_path=work_dir,
    )

    for region in regions_dict.keys():
        if region not in adata_tfbs.obsm:
            continue
        print(region, regions_dict[region])
        tfbs = pd.DataFrame(adata_tfbs.obsm[region], index=group_names)
        g = sns.clustermap(tfbs, cmap="RdBu_r", col_cluster=False)
        g.fig.suptitle(regions_dict[region])
        plt.savefig(os.path.join(plots_dir, f"TFBS_{regions_dict[region]}.png"), dpi=200)
        plt.close()

    # 5.8 Select pseudobulks to run LoRA on
    # For general datasets (not PBMC tutorial), just use ALL pseudobulk groups.
    selected_pbulks = list(group_names)
    print(f"[INFO] Using all {len(selected_pbulks)} pseudobulk groups for LoRA TFBS.")

    # Corresponding bigwig paths
    bigwig_dict = {
        sample: os.path.join(work_dir, f"{sample}_TFBS.bigwig")
        for sample in selected_pbulks
    }

    # Check if all pseudobulk TFBS bigwigs already exist
    all_tfbs_bigwigs_exist = (
        len(selected_pbulks) > 0
        and all(exists_nonempty(bigwig_dict[sample]) for sample in selected_pbulks)
    )

    if args.resume and all_tfbs_bigwigs_exist:
        print("[INFO] All pseudobulk TFBS bigwigs already exist; skip LoRA TFBS on regions_filt.")
    else:
        # 5.9 Run LoRA TFBS on regions_filt.bed for selected pseudobulks
        print("[INFO] Running LoRA TFBS on filtered regions (pseudobulks) ...")
        _ = scp.tl.seq_tfbs_seq2print(
            seq_attr_count=None,
            seq_attr_footprint=None,
            genome=printer.genome,
            region_path=regions_filt_bed,
            gpus=[0],
            model_type="lora",
            model_path=lora_model_path,
            lora_config=lora_cfg,          # reuse the already-patched config
            group_names=selected_pbulks,
            verbose=False,
            launch=True,
            return_adata=False,
            overwrite_seqattr=False,
            save_key="PBMC_LoRA",
            save_path=work_dir,
        )

    # 5.10 Motif scan + TFBS_scores from pseudobulk bigwigs
    print("[INFO] Scanning motifs on filtered regions ...")
    motif_sites_path = os.path.join(work_dir, "motif_sites.tsv")

    if args.resume and exists_nonempty(motif_sites_path):
        print(f"[INFO] motif_sites already exists, loading from {motif_sites_path}")
        motif_sites = pd.read_csv(motif_sites_path, sep="\t")
    else:
        print("[INFO] Scanning motifs on filtered regions ...")
        motifs = scp.motifs.FigR_Human_Motifs(genome=genome, bg=[0.25] * 4)
        motifs.prep_scanner()
        motif_sites = motifs.scan_motif(regions_filt, verbose=True, clean=True)
        motif_sites = pd.DataFrame(motif_sites)
        motif_sites.iloc[:, 2] = motif_sites.iloc[:, 1] + motif_sites.iloc[:, 8]
        motif_sites.iloc[:, 1] = motif_sites.iloc[:, 1] + motif_sites.iloc[:, 7]
        motif_sites = motif_sites.iloc[:, [0, 1, 2, 4]]
        motif_sites.columns = ["chrom", "start", "end", "TF"]
        motif_sites.to_csv(motif_sites_path, sep="\t", index=False)

    print("[INFO] Fetching TFBS bigwig scores for pseudobulks ...")
    TFBS_scores_path = os.path.join(work_dir, "TFBS_scores.csv")
    tfbs_scatter_path = os.path.join(work_dir, "TFBS_scatter.png")

    if args.resume and exists_nonempty(TFBS_scores_path) and exists_nonempty(tfbs_scatter_path):
        print(f"[INFO] TFBS_scores and TFBS_scatter already exist, skip bigwig fetching.")
    else:
        print("[INFO] Fetching TFBS bigwig scores for pseudobulks ...")

        TFBS_scores_list = []
        for sample in selected_pbulks:
            print(f"[INFO] Fetching TFBS scores for {sample} ...")
            scores = fetch_bw([motif_sites, bigwig_dict[sample], genome])
            TFBS_scores_list.append(scores)

        TFBS_scores = np.array(TFBS_scores_list).T
        TFBS_scores = pd.DataFrame(
            TFBS_scores,
            columns=[f"TFBS_{sample}" for sample in selected_pbulks],
        )

        TFBS_scores = pd.concat([motif_sites, TFBS_scores], axis=1)
        TFBS_scores_path = os.path.join(work_dir, "TFBS_scores.csv")
        TFBS_scores.to_csv(TFBS_scores_path, index=False)

        tfbs_scatter_path = os.path.join(work_dir, "TFBS_scatter.png")

        tfbs_cols = [c for c in TFBS_scores.columns if c.startswith("TFBS_")]

        if len(tfbs_cols) > 0:
            mean_tfbs = TFBS_scores[tfbs_cols].mean(axis=1)

            plt.figure(figsize=(6, 4))
            plt.scatter(
                np.arange(len(mean_tfbs)),
                mean_tfbs,
                s=2,
                alpha=0.5,
            )
            plt.xlabel("Motif site index")
            plt.ylabel("Mean TFBS score across pseudobulks")
            plt.title("TFBS scores per motif site (mean across pseudobulks)")
            plt.tight_layout()
            plt.savefig(tfbs_scatter_path, dpi=200)
            plt.close()
        else:
            plt.figure(figsize=(4, 3))
            plt.text(0.5, 0.5, "No TFBS_* columns found", ha="center", va="center")
            plt.axis("off")
            plt.tight_layout()
            plt.savefig(tfbs_scatter_path, dpi=200)
            plt.close()

    print("[INFO] Full single-cell PBMC-style pipeline (with LoRA) finished.")
    print(f"[INFO] Outputs under: {work_dir}")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()
