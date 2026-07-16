import os
import glob

script_dir = os.path.dirname(workflow.snakefile)
threads = int(config["threads"])
fdr = config["fdr"]
bin_size = config["bin_size"]
fitType = config["fitType"]
samplesheet = config["samples_sheet"]
samplesheet_tag = os.path.splitext(os.path.basename(samplesheet))[0]
outdir = f"differential_interactions_hicdcplus_{samplesheet_tag}"

def sample_paths_from_samplesheet(ss):
    paths = []
    with open(ss, "r") as f:
        next(f)  # skip header
        for line in f:
            if line.strip():
                paths.append(line.strip().split("\t")[0])
    return paths

def sample_conditions_from_samplesheet(ss):
    samples_comparison = []
    with open(ss, "r") as f:
        next(f)  # skip header
        for line in f:
            if line.strip():
                samples_comparison.append(line.strip().split("\t")[1])
    return samples_comparison
    
samples_comparison = sample_conditions_from_samplesheet(samplesheet)

rule all:
    input:
        f"{outdir}/hicdcdiff_sessionInfo.txt"

#############################################
# --- Differential analysis using hicdcplus Pipeline Rules ---
#############################################
rule hicdc_differential_interactions:
    input:
        samplesheet=samplesheet,
        sample_bedpes=lambda wildcards: sample_paths_from_samplesheet(samplesheet),
        # NEW: optional consensus file — set consensus_bedpe="" in config to skip
        consensus=config.get("consensus_bedpe", ""),
    output:
        pdf=f"{outdir}/hicdcdiff_sessionInfo.txt",
    params:
        outdir=outdir,
        samplesheet=samplesheet,
        fdr=fdr,
        binsize=bin_size,
        fitType=fitType,
        min_support=config.get("consensus_min_support", 2),
        # only pass --consensus flag when the file is actually provided
        consensus_arg=lambda wildcards, input:
            f"--consensus {input.consensus} --min_support {config.get('consensus_min_support', 2)}"
            if input.consensus and input.consensus != "" else "",
    log:
        f"logs/differential_interactions_hicdc.{samplesheet_tag}.log"
    benchmark:
        f"benchmarks/differential_interactions_hicdc.{samplesheet_tag}.txt"
    threads: 1
    conda:
        "envs/hicdc.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript {script_dir}/../scripts/hicdc_differential_interactions.R \
          --samplesheet {params.samplesheet} \
          --output_dir  {params.outdir}      \
          --fitType     {params.fitType}     \
          --binsize     {params.binsize}     \
          --fdr         {params.fdr}         \
          {params.consensus_arg}             \
          > {log} 2>&1
        """


#############################################
# --- Differential analysis using diffloop Pipeline Rules ---
#############################################
rule hichipper_diffloop_differential_interactions:
    input:
        control_files=glob.glob("hichipper/Hichipper.{samples_comparison[0]}*filt.intra.loop_counts.bedpe"),
        test_files=glob.glob("hichipper/Hichipper.{samples_comparison[1]}*filt.intra.loop_counts.bedpe")
    output:
        bedpe="differential_interactions_{samples_comparison[0]}_vs_{samples_comparison[1]}_diffloop/diffloop_pca_plot.tiff",
    params:
        indir="hichipper/",
        outdir="differential_interactions_{samples_comparison[0]}_vs_{samples_comparison[1]}_diffloop/",
        fdr=fdr,  
        control=samples_comparison[0],
        test=samples_comparison[1],
        control_files=" ".join([os.path.basename(f).replace(".filt.intra.loop_counts.bedpe", "") for f in glob.glob(f"hichipper/Hichipper.{samples_comparison[0]}*filt.intra.loop_counts.bedpe")]),
        test_files=" ".join([os.path.basename(f).replace(".filt.intra.loop_counts.bedpe", "") for f in glob.glob(f"hichipper/Hichipper.{samples_comparison[1]}*filt.intra.loop_counts.bedpe")]),
    log:
        "logs/hichipper_diffloop_differential_interactions.{samples_comparison[0]}_vs_{samples_comparison[1]}.log",
    benchmark:
        "benchmarks/hichipper_diffloop_differential_interactions.{samples_comparison[0]}_vs_{samples_comparison[1]}.txt",
    threads: 1
    conda:
        "envs/diffloop.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/hichipper.diffloop_differential_analysis.R \
          {params.control} {params.test} \
          "{params.control_files}" "{params.test_files}" \
          {params.indir} {params.outdir} \
          {params.fdr} > {log} 2>&1
        """

#############################################
# --- Differential analysis using Fithichip Pipeline Rules ---
#############################################
rule fithichip_differential_interactions:
    input:
        control_gz_files=glob.glob("hicdcplus/{samples_comparison[0]}*/Hicdcplus.significant_interactions.bedpe.gz"),
        test_gz_files=glob.glob("hicdcplus/{samples_comparison[1]}*/Hicdcplus.significant_interactions.bedpe.gz")
    output:
        bedpe="differential_interactions_{samples_comparison[0]}_vs_{samples_comparison[1]}_hicdcplus/dispersionplot.pdf",
    params:
        outdir="differential_interactions_{samples_comparison[0]}_vs_{samples_comparison[1]}_hicdcplus/",
        fdr=fdr,  
        control=samples_comparison[0],
        test=samples_comparison[1],
    log:
        "logs/fithichip_differential_interactions.{samples_comparison[0]}_vs_{samples_comparison[1]}.log"
    benchmark:
        "benchmarks/fithichip_differential_interactions.{samples_comparison[0]}_vs_{samples_comparison[1]}.txt"
    threads: 1
    conda:
        "envs/hicdcplus.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/DiffAnalysisHiChIP.R \
        —AllLoopList cat1_repl1_file,cat1_repl2_file,cat2_repl1_File,cat2_repl2_file \
        --ChrSizeFile TestData/chrom_hg19.sizes \
        —ChIPAlignFileList cat1_ChIPAlign.bam,cat2_ChIPAlign.bam \
        --FDRThr 0.01 \
        —BackgroundFDRThr 0.1 \
        --CovThr 25 \
        —OutDir /home/user/diffanalysis/outdir \
        --CategoryList CellLine1,CellLine2 \
        —ReplicaCount 2,2 \
        --ReplicaLabels1 R1,R2 \
        —ReplicaLabels2 R1,R2 \
        --FoldChangeThr 2 \
        —DiffFDRThr 0.05 \
        --bcv 0.4 > {log} 2>&1
        """



