import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"  # Limit to 2 GPUs (change if needed)
import torch
import pooch
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
import pyranges as pr
from tqdm.contrib.concurrent import *
from tqdm.auto import *
import anndata
import scanpy as sc
import json
import csv
import re
from sklearn.preprocessing import OneHotEncoder
import glob
import argparse
import snapatac2 as snap
import pyBigWig as pw
from pathlib import Path
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
def patch_scprinter_shift_detection_read_csv():
    try:
        import scprinter.shift_detection as _sd
        import pandas as pd
        import os
        from pathlib import Path

        if getattr(pd.read_csv, "_snakehichiptf_patched", False):
            print("[INFO] pandas.read_csv already patched")
            return

        _orig_read_csv = pd.read_csv

        def _patched_read_csv(filepath_or_buffer, *args, **kwargs):
            is_pathlike = isinstance(filepath_or_buffer, (str, os.PathLike, Path))
            if is_pathlike:
                p = os.fspath(filepath_or_buffer)
                print(f"[DEBUG] pd.read_csv called on: {p}")

                kwargs = dict(kwargs)

                if p.endswith((".gz", ".bgz")):
                    kwargs.setdefault("compression", "gzip")
                    print(f"[DEBUG] forcing gzip compression for: {p}")

            return _orig_read_csv(filepath_or_buffer, *args, **kwargs)

        _patched_read_csv._snakehichiptf_patched = True

        pd.read_csv = _patched_read_csv
        _sd.pd.read_csv = _patched_read_csv

        print("[INFO] Patched pandas.read_csv and scprinter.shift_detection.pd.read_csv")
    except Exception as e:
        print(f"[WARN] Failed to patch scprinter shift_detection read_csv: {e}")

def patch_scprinter_import_data_for_snapatac2():
    try:
        import scprinter.preprocessing as sp
        import snapatac2 as snap
        import inspect

        # Only patch when old API is missing
        if callable(getattr(snap.pp, "import_data", None)):
            print("[INFO] snapatac2.preprocessing.import_data exists; no scprinter patch needed")
            return

        # New SnapATAC2 API location
        import_fragments_fn = None
        if hasattr(snap.pp, "_import_data"):
            mod = snap.pp._import_data
            import_fragments_fn = getattr(mod, "import_fragments", None)

        if not callable(import_fragments_fn):
            raise RuntimeError("Could not find callable snapatac2.preprocessing._import_data.import_fragments")

        _orig_import_data = sp.import_data
        def _patched_import_data(*args, **kwargs):
            print("[WARN] Using patched scprinter.preprocessing.import_data -> snapatac2.import_fragments")
            print("[DEBUG] args =", args)
            print("[DEBUG] kwargs =", kwargs)
            print("[DEBUG] original kwargs keys =", list(kwargs.keys()))

            # old import_data positional args, inferred from your debug:
            # 0: fragment_file
            # 1: something currently None
            # 2: genome object
            # 3: some integer
            # 4: some integer
            # 5: output path
            fragment_file = args[0] if len(args) > 0 else kwargs.get("fragment_file")
            maybe_none    = args[1] if len(args) > 1 else None
            genome_obj    = args[2] if len(args) > 2 else kwargs.get("genome")
            arg3          = args[3] if len(args) > 3 else None
            arg4          = args[4] if len(args) > 4 else None
            out_path      = args[5] if len(args) > 5 else None

            print("[DEBUG] fragment_file =", repr(fragment_file))
            print("[DEBUG] maybe_none =", repr(maybe_none))
            print("[DEBUG] genome_obj =", repr(genome_obj))
            print("[DEBUG] genome_obj type =", type(genome_obj))
            print("[DEBUG] arg3 =", repr(arg3), "type =", type(arg3))
            print("[DEBUG] arg4 =", repr(arg4), "type =", type(arg4))
            print("[DEBUG] out_path =", repr(out_path))

            chrom_sizes = None

            # try to recover chrom_sizes from genome object
            if genome_obj is not None:
                for attr in ["chrom_sizes", "chromSize", "chromsizes", "_chrom_sizes"]:
                    if hasattr(genome_obj, attr):
                        chrom_sizes = getattr(genome_obj, attr)
                        print(f"[DEBUG] using genome_obj.{attr} for chrom_sizes")
                        break

            print("[DEBUG] chrom_sizes value =", repr(chrom_sizes))
            print("[DEBUG] chrom_sizes type =", type(chrom_sizes))
            if genome_obj is not None:
                print("[DEBUG] genome_obj dir sample =", [x for x in dir(genome_obj) if "chrom" in x.lower()])
            if chrom_sizes is None:
                raise ValueError(
                    "[PATCH DEBUG] chrom_sizes is still None. "
                    "Need to inspect genome_obj attributes and map the correct one."
                )

            mapped = {
                "fragment_file": fragment_file,
                "chrom_sizes": chrom_sizes,
                "min_num_fragments": kwargs.get("min_num_fragments", 1000),
                "sorted_by_barcode": kwargs.get("sorted_by_barcode", False),
            }

            print("[DEBUG] mapped kwargs keys =", list(mapped.keys()))
            print("[DEBUG] mapped =", mapped)

            return import_fragments_fn(**mapped)
        sp.import_data = _patched_import_data
        print("[WARN] Patched scprinter.preprocessing.import_data for SnapATAC2>=2.9")

    except Exception as e:
        print(f"[WARN] Failed to patch scprinter import_data: {e}")
        raise

# ------------------------------------------------------------------------

def auto_select_gpus():
    try:
        output = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL,
        ).decode().strip().splitlines()
        gpus = []
        for line in output:
            gpu_id, mem_used = line.split(",")
            gpus.append((int(mem_used.strip()), gpu_id.strip()))
        gpus.sort(key=lambda x: x[0])
        if len(gpus) <= 2:
            return [g[1] for g in gpus]
        return [gpus[0][1], gpus[1][1]]
    except Exception:
        return ["0", "1"]  # fallback

# Patch SnapATAC2 TSSe
def _dummy_tsse(adata, gene_anno):
    return np.zeros(adata.n_obs, dtype=float)
snap.metrics.tsse = _dummy_tsse

def exists_nonempty(path: str) -> bool:
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0

# ====================== ROBUST fetch_bw (fixes "Invalid interval bounds!") ======================
def fetch_bw(args):
    import pyBigWig as pw
    import numpy as np
    from tqdm import tqdm
    
    TFBS, bw_path, genome = args
    
    # Motif sites from scPrinter (always use "chr" prefix)
    motif_chroms = np.array(TFBS["chrom"])
    starts = np.array(TFBS["start"])
    ends = np.array(TFBS["end"])
    
    res_all = {}
    with pw.open(bw_path) as f:
        bw_chroms = f.chroms()  # Get actual chromosome names in this bigWig
        
        if bw_chroms is None:
            print(f"[WARN] BigWig {bw_path} has no chromosomes — returning NaN for all sites")
            return [np.nan] * len(TFBS)
        
        # Detect if this bigWig uses "chr" prefix
        has_chr = any(k.startswith("chr") for k in bw_chroms.keys())
        
        # Pre-load full chromosome signals
        for bw_chrom, length in tqdm(bw_chroms.items(), desc=f"Loading chroms from {os.path.basename(bw_path)}", leave=False):
            # Skip non-standard or unwanted contigs
            if bw_chrom in ["chrY", "Y"] or bw_chrom.startswith(("chrUn_", "random", "EBV", "chrEBV", "alt", "fix")):
                continue
            try:
                signal = f.values(bw_chrom, 0, length, numpy=True)
                res_all[bw_chrom] = signal
            except RuntimeError as e:
                print(f"[WARN] Failed to read {bw_chrom} from {bw_path}: {e}")
                continue
    
    # Extract mean signal for each motif site
    vs = []
    for chrom, start, end in zip(tqdm(motif_chroms, desc="Extracting TFBS signals", mininterval=1), starts, ends):
        # Primary attempt: direct match
        key = chrom
        if key not in res_all:
            # Fallback: toggle "chr" prefix based on bigWig style
            if has_chr and not chrom.startswith("chr"):
                key = "chr" + chrom
            elif not has_chr and chrom.startswith("chr"):
                key = chrom[3:]  # remove "chr"
        
        if key in res_all:
            segment = res_all[key][start:end]
            vs.append(np.nanmean(segment))
        else:
            vs.append(np.nan)  # Motif on scaffold/unmapped chrom
    
    return vs
# ===============================================================================================

def select_regions_with_footprints(peaks_bed, footprints_bed, out_bed, n_regions=10):
    peaks = pd.read_csv(peaks_bed, sep="\t", header=None, usecols=[0,1,2], names=["chrom","start","end"])
    footprints = pd.read_csv(footprints_bed, sep="\t", header=None, usecols=[0,1,2], names=["chrom","start","end"])
    peaks_pr = pr.PyRanges(peaks)
    fp_pr = pr.PyRanges(footprints)
    overlapped = peaks_pr.join(fp_pr)
    if len(overlapped) == 0:
        raise ValueError(f"No overlap between {peaks_bed} and {footprints_bed}.")
    df = overlapped.df[["Chromosome", "Start", "End"]].drop_duplicates()
    df = df.sort_values(["Chromosome", "Start", "End"])
    selected = df.head(n_regions)
    selected.to_csv(out_bed, sep="\t", header=False, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main_dir", required=True)
    parser.add_argument("--entity", required=True)
    parser.add_argument("--project", default="scPrinter_ATAC")
    parser.add_argument("--regions_bed", default=None)
    parser.add_argument("--footprints_bed", default=None)
    parser.add_argument("--n_regions", type=int, default=10)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    main_dir = args.main_dir
    entity = args.entity
    project = args.project
    regions_bed_arg = args.regions_bed
    footprints_bed = args.footprints_bed
    n_regions = args.n_regions

    work_dir = os.path.join(main_dir, "seq2print")
    os.makedirs(work_dir, exist_ok=True)

    frag_dir = f"{main_dir}/fragments"
    os.makedirs(frag_dir, exist_ok=True)
    frag_files = sorted([
        os.path.join(frag_dir, f)
        for f in os.listdir(frag_dir)
        if f.endswith("_fragments.tsv.gz")
    ])
    #frag_files = sorted([os.path.join(frag_dir, f) for f in os.listdir(frag_dir) if re.search("fragments.tsv.gz", f)])
    samples = [os.path.basename(f).replace("_fragments.tsv.gz", "") for f in frag_files]

    gpu_ids = auto_select_gpus()
    print(f"[INFO] Using GPUs: {gpu_ids}")
    print(f"[DEBUG] main_dir = {main_dir}")
    print(f"[DEBUG] frag_dir = {frag_dir}")
    print(f"[DEBUG] frag_files = {frag_files}")
    print(f"[DEBUG] samples = {samples}")

    patch_scprinter_shift_detection_read_csv()
    patch_scprinter_import_data_for_snapatac2()
    # Import fragments
    printer = scp.pp.import_fragments(
        path_to_frags=frag_files,
        barcodes=[None] * len(frag_files),
        savename=os.path.join(work_dir, "ATAC_scprinter.h5ad"),
        genome=scp.genome.hg38,
        min_num_fragments=1000,
        min_tsse=7,
        sorted_by_barcode=False,
        low_memory=False,
    )

    # Call peaks
    cleaned_peaks_path = os.path.join(work_dir, "seq2print_cleaned_narrowPeak.bed")
    if args.resume and exists_nonempty(cleaned_peaks_path):
        print(f"[INFO] Resume: reusing cleaned peaks bed")
    else:
        scp.pp.call_peaks(printer=printer, frag_file=frag_files, cell_grouping=[None], group_names=["all"], preset="seq2PRINT", overwrite=False)
        pd.DataFrame(printer.uns["peak_calling"]["all_cleaned"][:]).to_csv(cleaned_peaks_path, sep="\t", header=False, index=False)

    chromvar_regions_path = os.path.join(work_dir, "regions.bed")
    if args.resume and exists_nonempty(chromvar_regions_path):
        print(f"[INFO] Resume: reusing chromvar regions bed")
    else:
        scp.pp.call_peaks(printer=printer, frag_file=frag_files, cell_grouping=[None], group_names=["chromvar_all"], preset="chromvar", overwrite=False)
        pd.DataFrame(printer.uns["peak_calling"]["chromvar_all_cleaned"][:]).to_csv(chromvar_regions_path, sep="\t", header=False, index=False)

    # Peak matrix and filtering
    adata = scp.pp.make_peak_matrix(printer, regions=chromvar_regions_path, region_width=300, sparse=True)
    adata.write(os.path.join(work_dir, "cell_peak.h5ad"))

    regions = pd.read_csv(chromvar_regions_path, sep="\t", header=None)
    adata = anndata.read_h5ad(os.path.join(work_dir, "cell_peak.h5ad"))
    peak_depth = np.array(np.sum(adata.X, axis=0)).squeeze()
    regions_filt = regions.iloc[np.where(peak_depth > 200)[0], :]

    # OOM protection: limit number of regions
    MAX_REGIONS = 30000  # Adjust this value if needed
    if len(regions_filt) > MAX_REGIONS:
        print(f"[INFO] Limiting regions from {len(regions_filt)} to {MAX_REGIONS} to prevent GPU OOM")
        regions_filt = regions_filt.head(MAX_REGIONS)

    regions_filt_path = os.path.join(work_dir, "regions_filt.bed")
    regions_filt.to_csv(regions_filt_path, sep="\t", header=False, index=False)

    # Train seq2PRINT models
    bulk_barcodes = np.array(printer.obs_names)
    configs_dir = os.path.join(work_dir, "configs")
    os.makedirs(configs_dir, exist_ok=True)
    region_filename_for_model = os.path.basename(cleaned_peaks_path)
    os.makedirs(os.path.join(work_dir, "temp"), exist_ok=True)
    os.makedirs(os.path.join(work_dir, "model"), exist_ok=True)

    model_path_dict = {}
    failed_samples = []
    train_items = []

    for i, sample in enumerate(samples):
        existing_pattern = os.path.join(work_dir, "model", f"ATAC_{sample}*.pt")
        existing_models = glob.glob(existing_pattern)
        if args.resume and existing_models:
            print(f"[INFO] Found existing model for {sample}, skipping training")
            model_path_dict[sample] = existing_models[0]
            continue
        train_items.append((i, sample))

    if train_items:
        max_workers = min(len(gpu_ids), len(train_items))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for i, sample in train_items:
                gpu = gpu_ids[i % len(gpu_ids)]
                cfg_path = os.path.join(configs_dir, f"ATAC_{sample}_fold0.JSON")
                if os.path.exists(cfg_path):
                    os.remove(cfg_path)
                scp.tl.seq_model_config(
                    printer=printer,
                    region_path=region_filename_for_model,
                    cell_grouping=[bulk_barcodes],
                    group_names=[sample],
                    genome=printer.genome,
                    fold=0,
                    overwrite_bigwig=False,
                    model_name="ATAC_" + sample,
                    additional_config={
                        "notes": "NPC bulk per-sample",
                        "tags": ["ATAC", sample],
                        "entity": entity,
                        "project": project,
                    },
                    config_save_path=cfg_path,
                )
                # ===== PATCH signals PATH =====
                with open(cfg_path) as f:
                    cfg = json.load(f)
                sig = cfg.get("signals")
                prefix = work_dir.rstrip("/") + "/"
                if isinstance(sig, str) and sig.startswith(prefix):
                    cfg["signals"] = sig[len(prefix):]
                with open(cfg_path, "w") as f:
                    json.dump(cfg, f, indent=2)
                sample_temp_dir = os.path.join(work_dir, "temp", f"train_{sample}")
                os.makedirs(sample_temp_dir, exist_ok=True)
                fut = executor.submit(
                    scp.tl.launch_seq2print,
                    model_config_path=cfg_path,
                    temp_dir=sample_temp_dir,
                    model_dir=os.path.join(work_dir, "model"),
                    data_dir=work_dir,
                    gpus=gpu,
                    verbose=True,
                    launch=True,
                )
                futures[fut] = sample
            for fut in as_completed(futures):
                sample = futures[fut]
                try:
                    fut.result()
                except Exception as e:
                    print(f"[WARN] Training failed for {sample}: {e}")
                    failed_samples.append(sample)
                    continue
                pattern = os.path.join(work_dir, "model", f"ATAC_{sample}*.pt")
                candidates = glob.glob(pattern)
                if candidates:
                    model_path_dict[sample] = candidates[0]
                else:
                    failed_samples.append(sample)

    if not model_path_dict:
        print("[ERROR] No models trained — writing dummy outputs")
        TFBS_scores_path = os.path.join(work_dir, "TFBS_scores.csv")
        TFBS_scatter_path = os.path.join(work_dir, "TFBS_scatter.png")
        pd.DataFrame(columns=["chrom","start","end","TF"]).to_csv(TFBS_scores_path, index=False)
        plt.figure(figsize=(4,3))
        plt.text(0.5, 0.5, "No model trained", ha="center", va="center")
        plt.axis("off")
        plt.savefig(TFBS_scatter_path, dpi=150)
        plt.close()
        printer.close()
        return

    # Region selection for visualization
    if regions_bed_arg and regions_bed_arg.strip():
        regions_test_path = regions_bed_arg
    else:
        regions_test_path = os.path.join(work_dir, "regions_test.bed")
        if footprints_bed and footprints_bed.strip():
            select_regions_with_footprints(cleaned_peaks_path, footprints_bed, regions_test_path, n_regions)
        else:
            pd.read_csv(cleaned_peaks_path, sep="\t", header=None).head(n_regions).to_csv(regions_test_path, sep="\t", header=False, index=False)

    # TFBS inference
    torch.cuda.empty_cache()

    def get_tfbs_bigwig_path(sample):
        return os.path.join(work_dir, f"tfbs_{sample}", f"{sample}_TFBS.bigwig")

    tfbs_bw_dict = {s: get_tfbs_bigwig_path(s) for s in model_path_dict}

    items_to_process = [(s, p) for s, p in tfbs_bw_dict.items() if not (args.resume and exists_nonempty(p))]

    if items_to_process:
        max_workers = min(2, len(items_to_process))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for i, (sample, bw_path) in enumerate(items_to_process):
                gpu = str(i % 2)
                sample_work_dir = os.path.join(work_dir, f"tfbs_{sample}")
                os.makedirs(sample_work_dir, exist_ok=True)
                sample_dummy_tfbs = os.path.join(tfbs_dir, f"{sample}_TFBS.bigwig")
                #sample_dummy_tfbs = os.path.join(sample_work_dir, "dummy_TFBS.bigwig")
                lora_config = json.load(open(os.path.join(configs_dir, f"ATAC_{sample}_fold0.JSON")))
                fut = executor.submit(
                    scp.tl.seq_tfbs_seq2print,
                    seq_attr_count=None,
                    seq_attr_footprint=None,
                    genome=printer.genome,
                    region_path=regions_filt_path,
                    gpus=gpu,
                    model_type="seq2print",
                    model_path=model_path_dict[sample],
                    lora_config=lora_config,
                    group_names=[sample],
                    verbose=True,
                    launch=True,
                    return_adata=False,
                    overwrite_seqattr=True,
                    post_normalize=True,
                    save_key=sample,
                    save_path=sample_work_dir,
                )
                futures[fut] = (sample, bw_path, sample_dummy_tfbs)
            for fut in as_completed(futures):
                sample, bw_path, sample_dummy_tfbs = futures[fut]
                try:
                    fut.result()
                except Exception as e:
                    print(f"[ERROR] TFBS inference failed for {sample}: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
                # if not exists_nonempty(sample_dummy_tfbs):
                #     raise FileNotFoundError(f"No dummy_TFBS.bigwig for {sample}")
                shutil.copy2(sample_dummy_tfbs, bw_path)

    # Motif scanning and scores
    TFBS_scores_path = os.path.join(work_dir, "TFBS_scores.csv")
    if args.resume and exists_nonempty(TFBS_scores_path):
        print("[INFO] Loading existing TFBS_scores.csv")
        TFBS_scores = pd.read_csv(TFBS_scores_path)
    else:
        print("[INFO] Running motif scanning and TFBS signal extraction")
        motifs = scp.motifs.FigR_Human_Motifs(genome=printer.genome, bg=[0.25]*4)
        motifs.prep_scanner()
        motif_sites = pd.DataFrame(motifs.scan_motif(regions_filt, verbose=True, clean=True))
        motif_sites.iloc[:, 2] = motif_sites.iloc[:, 1] + motif_sites.iloc[:, 8]
        motif_sites.iloc[:, 1] = motif_sites.iloc[:, 1] + motif_sites.iloc[:, 7]
        motif_sites = motif_sites.iloc[:, [0,1,2,4]]
        motif_sites.columns = ["chrom","start","end","TF"]
        
        # Prepare arguments: [motif_sites_df, bigwig_path, genome]
        args_bw = [[motif_sites, tfbs_bw_dict[s], printer.genome] for s in model_path_dict]
        
        import multiprocessing as mp
        with mp.Pool(4) as pool:
            scores_list = list(pool.imap(fetch_bw, args_bw))
        
        scores_df = pd.DataFrame(np.array(scores_list).T, columns=[f"TFBS_{s}" for s in model_path_dict])
        TFBS_scores = pd.concat([motif_sites.reset_index(drop=True), scores_df], axis=1)
        TFBS_scores.to_csv(TFBS_scores_path, index=False)
        print(f"[INFO] Saved TFBS_scores.csv with {len(TFBS_scores)} motif instances")

    # Scatter plots
    TFBS_scatter_path = os.path.join(work_dir, "TFBS_scatter.png")
    if not (args.resume and exists_nonempty(TFBS_scatter_path)):
        if len(model_path_dict) >= 2:
            s1, s2 = list(model_path_dict.keys())[:2]
            fig, ax = plt.subplots(2,2, figsize=(7,6))
            for idx, tf in enumerate(["CTCF", "SPI1", "KLF1", "RUNX3"]):
                sub = TFBS_scores[TFBS_scores.TF == tf]
                r, c = divmod(idx, 2)
                ax[r][c].scatter(sub[f"TFBS_{s1}"], sub[f"TFBS_{s2}"], s=0.01)
                ax[r][c].set_xlabel(s1)
                ax[r][c].set_ylabel(s2)
                ax[r][c].set_title(tf)
            plt.tight_layout()
            plt.savefig(TFBS_scatter_path, dpi=300)
            plt.close()
        else:
            plt.figure(figsize=(4,3))
            plt.text(0.5, 0.5, "Only one sample", ha="center", va="center")
            plt.axis("off")
            plt.savefig(TFBS_scatter_path, dpi=150)
            plt.close()

    printer.close()

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()
