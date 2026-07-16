import os
import yaml
import glob

#############################################
# Configuration
#############################################
script_dir = os.path.dirname(workflow.snakefile)

# Load configuration
pipeline = config["pipeline"]
input_dir = config["input_dir"]
organism = config["organism"]
genome_dir = config["genome_dir"]
downsample_size = config["downsample_size"]
samples = config["samples"]
reads = config["reads"]
hicpro_params = config["hicpro_params"]
fithichip_BiasType = config["fithichip_BiasType"]
hichipper_params = config["hichipper_params"]
hicdc_params = config["hicdc_params"]
restriction_enzyme = config["restriction_enzyme"]
bin_size = config["bin_size"]
binning_range = config["binning_range"]
generate_hic = config["generate_hic"]
mapq = config["mapq"]
length_cutoff = config["length_cutoff"]
threads = int(config["threads"] / 2)
maps_model = config["maps_model"]
maps_sex_chroms = config["maps_sex_chroms"]
macs2_peaks = config["macs2_peaks"]
optical_duplicate_distance = config["optical_duplicate_distance"]
maps_count_cutoff = config["maps_count_cutoff"]
maps_ratio_cutoff = config["maps_ratio_cutoff"]
fdr = config["fdr"]
resolution = bin_size // 1000

# Parse sample comparison from command line (e.g., --samples_comparison "SampleA SampleB")
samples_comparison = config["samples_comparison"]
samples_comparison = samples_comparison.split()

#############################################
# Final target
#############################################
# Define pipelines and their associated outputs
PIPELINE_OUTPUTS = {
    "Maps": [
        "maps_feather/{sample}/{sample}.long.intra.bedpe",
        "maps/{sample}/Maps.{sample}." + str(resolution) + "k.sig3Dinteractions.bedpe"
    ],
    "Fithichip": [
        "fithichip/{sample}/FitHiChIP_alltoall_BinSize" + str(bin_size) + "/Plots/FitHiChIP.{sample}.interactions.alltoall.significant.Dist_CC.png",
        "fithichip/{sample}/FitHiChIP_peaktoall_BinSize" + str(bin_size) + "/Plots/FitHiChIP.{sample}.interactions.peaktoall.significant.merged.Dist_CC.png"
    ],
    "Hichipper": [
        "hichipper/Hichipper.{sample}.significant.interactions.bedpe",
        "hichipper/Hichipper.{sample}_qc_report.html"
    ],
    "Hicdcplus": [
        "hicdcplus/{sample}/Hicdcplus.{sample}.significant_interactions.bedpe"
    ]
    "Consensus": [                                             
        "consensus/{sample}.consensus_interactions.bedpe"
    ]

}
rule all:
    input:
        "config_used.yaml",
        "qc/summary_hichip_valid_pairs_barplot.pdf",
        "qc/summary_hichip_frip_barplot.pdf",
        "quantiles/quantile_dot_plot_hichip.pdf",
        "quantiles/quantile_dot_plot_hichip_over1.pdf",
        # Prepared fastqs
        expand("processedFastqs/{sample}/{sample}_{read}.fastq.gz", read=reads, sample=samples),
        # Hicpro mapping
        expand("hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs", sample=samples), #Hicpro result
        # MACS2 peak calling
        expand("peaks/{sample}_peaks.narrowPeak", sample=samples),
        # Pipeline-specific outputs
        lambda wildcards: (
            # If pipeline is "All", include all pipeline outputs
            sum([expand(pattern, sample=samples) for pattern in sum(PIPELINE_OUTPUTS.values(), [])], [])
            if config.get("pipeline", "All") == "All"
            # Otherwise, include only the outputs for the specified pipeline
            else expand(PIPELINE_OUTPUTS[config["pipeline"]], sample=samples)
        ),
        # Quality control
        "qc/multiqc_report.html", 
        # Differential analysis output if samples_comparison is provided
        #lambda wc: expand("hicdc_diff_output/{sampleA}_vs_{sampleB}_diff_interactions.bedpe",
        #                  sampleA=samples_comparison[0], sampleB=samples_comparison[1]) if samples_comparison else []


onstart:
    shell("mkdir -p originalFastqs processedFastqs fastqc maps_feather maps logs benchmarks hichipper")

#############################################
# Config file output
#############################################
rule show_config:
    output: "config_used.yaml"
    run:
        import yaml
        with open(output[0], "w") as f:
            yaml.dump(config, f, sort_keys=True)

#############################################
# --- FASTQ Preprocessing Rules ---
#############################################
rule link_fastqs:
    input:
        fastqs=lambda wc: os.path.join(input_dir, f"{wc.sample}_{wc.read}.fastq.gz")
    output:
        "originalFastqs/{sample}_{read}.fastq.gz"
    log:
        "logs/link_fastqs.{sample}_{read}.log"
    benchmark:
        "benchmarks/link_fastqs.{sample}_{read}.txt"
    threads: threads
    conda:
        "envs/basic.yaml"
    shell:
        """
        ln -s {input.fastqs} {output} 2> {log}
        """

rule downsample:
    input:
        "originalFastqs/{sample}_{read}.fastq.gz"
    output:
        "downsampledFastqs/{sample}_{read}.fastq.gz"
    params:
        size=downsample_size
    log:
        "logs/downsample.{sample}_{read}.log"
    benchmark:
        "benchmarks/downsample.{sample}_{read}.txt"
    threads: threads
    conda:
        "envs/basic.yaml"
    shell:
        """
        seqtk sample -s100 {input} {params.size} | gzip > {output} 2> {log}
        """

# HELPER FUNCTION TO CHOOSE FASTQ SOURCE FOR HiCPro
def get_fastqs(wildcards):
    if int(downsample_size) > 0:
        return expand("downsampledFastqs/{sample}_{read}.fastq.gz", sample=wildcards.sample, read=reads)
    else:
        return expand("originalFastqs/{sample}_{read}.fastq.gz", sample=wildcards.sample, read=reads)

rule prepare_fastqs:
    input:
        #fastqs=lambda wildcards: get_fastqs(wildcards),
        fastqs=expand("originalFastqs/{{sample}}_{read}.fastq.gz", read=reads)
    output:
        expand("processedFastqs/{{sample}}/{{sample}}_{read}.fastq.gz", read=reads)
    params:
        raw_dir="processedFastqs/{sample}"
    log:
        "logs/prepare_fastqs.{sample}.log"  # Simplified log naming
    benchmark:
        "benchmarks/prepare_fastqs.{sample}.txt"
    shell:
        """
        mkdir -p {params.raw_dir} 2>> {log}
        for fq in {input.fastqs}; do
            ln -s ../../$fq {params.raw_dir}/$(basename $fq) 2>> {log}
        done
        """

###############################################################################
# --- HiC-Pro v3.1.0 for mapping and HiC ICE matrix generation ---
###############################################################################

# Rule to run HiC-Pro
rule hicpro_run:
    input:
        fastqs=expand("processedFastqs/{{sample}}/{{sample}}_{read}.fastq.gz", read=reads),
    output:
        bam="hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.bam",
        valid_pairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
        ice_matrix="hicpro/{sample}/hic_results/matrix/{sample}/iced/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_iced.matrix",
        matrix="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+".matrix",
        bin_intervals = "hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_abs.bed",
    params:
        config= genome_dir + "/organisms/"+organism+"/"+restriction_enzyme+"/"+str(bin_size)+"/"+organism+".txt" if hicpro_params == "NA" else hicdc_params,
        original_pairs = "hicpro/{sample}/hic_results/data/{sample}/{sample}_original.allValidPairs",
        out_dir = "hicpro/{sample}/hic_results/matrix/{sample}",
        input_dir = "hicpro/{sample}/hic_results/data/{sample}",
        raw_matrix = "hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size),
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
    log:
        "logs/hicpro_run.{sample}.log"
    benchmark:
        "benchmarks/hicpro_run.{sample}.txt"
    threads: threads
    conda:
        "envs/hicpro.yaml"
    shell:
        """
        {script_dir}/../scripts/hicpro-3.1.0/bin/HiC-Pro -i processedFastqs/{wildcards.sample} -o hicpro/{wildcards.sample} -c {params.config} > {log} 2>&1
        # Downsample if specified
        if [ {downsample_size} -gt 0 ]; then
            echo "Downsampling" >> {log}
            mv {output.valid_pairs} {params.original_pairs}  2>> {log}
            echo "Step1" >> {log}
            shuf --random-source=<(yes 1000) -n {downsample_size} {params.original_pairs} > {output.valid_pairs} 2>> {log}
            
            mv {params.out_dir} {params.out_dir}_original 2>> {log}
            mkdir -p {params.out_dir} 2>> {log}
            mkdir -p {params.out_dir}/raw/{bin_size} 2>> {log}
            mkdir -p {params.out_dir}/iced/{bin_size} 2>> {log}
            echo "Step2" >> {log}
            cat {output.valid_pairs} | {script_dir}/../scripts/hicpro-3.1.0/scripts/build_matrix --matrix-format upper --binsize {bin_size} --chrsizes {params.chromosome_sizes} --ifile /dev/stdin --oprefix {params.raw_matrix}  >> {log} 2>&1
            echo "Step3" >> {log}
            python {script_dir}/../scripts/hicpro-3.1.0/scripts/ice {output.matrix} --results_filename {output.ice_matrix}  --output-bias 1 --eps 0.1  --max_iter 100 --filter_high_counts_perc 0 --filter_low_counts_perc 0.02 --remove-all-zeros-loci >> {log} 2>&1
        else
            echo "No downsampling" >> {log}
        fi
        """

###############################################################################
# --- MACS2 for peak calling ---
###############################################################################

rule macs2_callpeak:
    input:
        bam="hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.bam",
    output:
        peaks="peaks/{sample}_peaks.narrowPeak",
    params:
        outdir="peaks",
        macs2_string="-q 0.1 --extsize 200 --nomodel --keep-dup all",  
        genome = lambda wildcards: "hs" if organism.startswith("hg") else organism[:2],
        name="{sample}",
    log:
        "logs/macs2_callpeak.{sample}.log",
    benchmark:
        "benchmarks/macs2_callpeak.{sample}.txt",
    threads: 1
    conda:
        "envs/hichipper.yaml"
    shell:
        """
        macs2 callpeak -t {input.bam} -f BAM -g {params.genome} {params.macs2_string} --outdir {params.outdir} --name {params.name} > {log} 2>&1
        """

#############################################
# --- Hicdc+ Pipeline Rules ---
#############################################
# Rule to identify significant interactions
rule hicdc_significant_interactions:
    input:
        valid_pairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
    output:
        bedpe="hicdcplus/{sample}/Hicdcplus.{sample}.significant_interactions.bedpe",
        hic="hicdcplus/{sample}/Hicdcplus.{sample}.all_interactions.hic"
    params:
        bintolen=expand(genome_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.hicdcplus.bintolen.txt.gz", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size),
        outdir="hicdcplus/{sample}",
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
        hicdc_params="" if hicdc_params == "NA" else hicdc_params,
        sample="{sample}",
        bedpe_gz="hicdcplus/{sample}/Hicdcplus.{sample}.significant_interactions.bedpe.gz",
    log:
        "logs/hicdc_significant_interactions.{sample}.log"
    benchmark:
        "benchmarks/hicdc_significant_interactions.{sample}.txt"
    threads: threads
    conda:
        "envs/hicdc.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/hicdc_significant_interactions.R \
          -i {input.valid_pairs} \
          -o {params.outdir} \
          -b {params.bintolen} \
          -c {params.chromosome_sizes} \
          -f {fdr} -s {params.sample} > {log} 2>&1
        gunzip -c {params.bedpe_gz} \
          | awk 'NR==1 {{print "#" $0; next}} {{print}}' \
          > {output.bedpe}
        """

#############################################
# --- Hichipper Pipeline Rules ---
#############################################
rule hichipper_run:
    input:
        allvalidpairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
        peaks="peaks/{sample}_peaks.narrowPeak",
    output:
        loops="hichipper/Hichipper.{sample}.filt.intra.loop_counts.bedpe",
    params:
        allvalidpairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
        peaks="peaks/{sample}_peaks.narrowPeak" if macs2_peaks == "NA" else macs2_peaks,
        outdir="{sample}",
        sample="{sample}",
        res=expand(genome_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_resolution_fragments.bed", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size),
        min_dist=5000,          # Minimum distance for interactions (adjust as needed)
        max_dist=binning_range,       # Maximum distance for interactions (adjust as needed)
        hichipper_params = "" if hichipper_params == "NA" else hichipper_params,
    log:
        "logs/hichipper_run.{sample}.log"
    benchmark:
        "benchmarks/hichipper_run.{sample}.txt"
    threads: threads
    conda:
        "envs/hichipper.yaml"
    shell:
        """
        if [ -d "{params.outdir}" ]; then
            rm -rf {params.outdir}
        fi 
        hichipper --input-vi {params.allvalidpairs} --peaks {params.peaks} --out {params.outdir} --restriction-frags {params.res} --min-dist {params.min_dist} --max-dist {params.max_dist} --skip-diffloop --basic-qc {params.hichipper_params} call > {log} 2>&1
        mv {params.outdir}/one.anchors.bed {params.outdir}/Hichipper.{params.sample}.anchors.bed
        mv {params.outdir}/one.filt.intra.loop_counts.bedpe {params.outdir}/Hichipper.{params.sample}.filt.intra.loop_counts.bedpe
        mv {params.outdir}/one.inter.loop_counts.bedpe {params.outdir}/Hichipper.{params.sample}.inter.loop_counts.bedpe
        mv {params.outdir}/one.intra.loop_counts.bedpe {params.outdir}/Hichipper.{params.sample}.intra.loop_counts.bedpe
        mv {params.outdir}/one.stat {params.outdir}/Hichipper.{params.sample}.stat
        mkdir -p hichipper
        mv {params.outdir}/* hichipper/
        rm -rf {params.outdir}
        """

rule hichipper_qc_report:
    input:
        loops="hichipper/Hichipper.{sample}.filt.intra.loop_counts.bedpe",
    output:
        qc_report="hichipper/Hichipper.{sample}_qc_report.html",
        rds="hichipper/Hichipper.{sample}.significant.interactions.rds",
        bedpe="hichipper/Hichipper.{sample}.significant.interactions.bedpe",     
    params:
        outputfile=os.getcwd()+"/hichipper/Hichipper.{sample}_qc_report.html",
        rmd=script_dir+"/../scripts/qcReport_make.Rmd",  # Path to your RMD file
        outdir="hichipper",
        scriptdir=script_dir,
        sample="Hichipper.{sample}",
        cwd=os.getcwd(),  # Dynamically gets the current working directory
        version="0.7.7",  # Adjust based on your HiChIPper version
    log:
        "logs/hichipper_qc_report.{sample}.log",
    benchmark:
        "benchmarks/hichipper_qc_report.{sample}.txt",
    threads: 1
    conda:
        "envs/diffloop.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/hichipper.diffloop_work.R {fdr} {maps_count_cutoff} {params.outdir} {params.sample} > {log} 2>&1        
        Rscript -e "rmarkdown::render('{params.rmd}', output_file = '{params.outputfile}', params = list(outdir = '{params.outdir}', scriptdir = '{params.scriptdir}', samples = '{params.sample}', cwd = '{params.cwd}', version = '{params.version}'))" >> {log} 2>&1
        """

#############################################
# --- FitHiChip Pipeline Rules ---
#############################################
rule fithichip_initial_interactions:
    input:
        valid_pairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
        matrix="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+".matrix",
    output:
        interactions="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.initial.bed",
    params:
        outdir="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"",
        bin_interval="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_abs.bed",
    log:
        "logs/fithichip_initial_interactions.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/InteractionHicPro.r {params.bin_interval} {input.matrix} {output.interactions} > {log} 2>&1
        """

rule fithichip_filter_interactions:
    input:
        interactions="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.initial.bed",
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
    output:
        filtered="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.cis.interactions.DistThr.bed",
    params:
        outdir="fithichip/{sample}/Matrix_BinSize"+str(bin_size),
        low_dist=5000,
        upp_dist=binning_range,
        circular=0,
    log:
        "logs/fithichip_filter_interactions.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/Filt_Loop_DistThr.r {input.interactions} {output.filtered} \
            {params.low_dist} {params.upp_dist} {params.circular} {input.chromosome_sizes} > {log} 2>&1
        """

rule fithichip_coverage:
    input:
        interactions="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.initial.bed",
        peaks="peaks/{sample}_peaks.narrowPeak",
    output:
        coverage="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.coverage.bed",
    params:
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
        interactions="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.initial.bed",
        peaks="peaks/{sample}_peaks.narrowPeak" if macs2_peaks == "NA" else macs2_peaks,
        outdir="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"",
    log:
        "logs/fithichip_coverage.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {script_dir}/../scripts/fithichip/CoverageBin.r --InpFile {input.interactions} --PeakFile {params.peaks} \
            --BinSize {bin_size} --ChrSizeFile {params.chromosome_sizes} --OutFile {output.coverage} > {log} 2>&1
        """

rule fithichip_bias_coverage:
    input:
        coverage="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.coverage.bed",
    output:
        coverage_bias="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.coverage.bed",
    params:
        features=expand(genome_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.fithichip.bed", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size),
        features_sorted="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/restriction_enzyme.fithichip.bed.sorted",
    log:
        "logs/fithichip_bias_coverage.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/BiasCalc.r --CoverageFile {input.coverage} --OutFile {output.coverage_bias}  > {log} 2>&1
        #Rscript {script_dir}/../scripts/fithichip/BiasCalc.r --CoverageFile {input.coverage} --OutFile {output.coverage_bias}.temp  > {log} 2>&1
        #awk 'NR>1' {output.coverage_bias}.temp | sort -k1,1 -k2,2n > {output.coverage_bias}.sorted
        #sort -k1,1 -k2,2n {params.features} > {params.features_sorted}
        #bedtools map -c 4,5 -o mean -null '0' -a {output.coverage_bias}.sorted -b {params.features_sorted} > {output.coverage_bias}.temp1
        #bedtools map -a {output.coverage_bias}.temp1 -b {params.features_sorted} -c 4 -o count -null '0' > {output.coverage_bias}
        #rm {output.coverage_bias}.temp {output.coverage_bias}.sorted {output.coverage_bias}.temp1
        """

rule fithichip_bias_ice:
    input:
        coverage="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.coverage.bed",
    output:
        coverage_bias="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.ice.bed",
    params:
        features=expand(genome_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.fithichip.bed", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size),
        bin_intervals = "hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_abs.bed",
        ice="hicpro/{sample}/hic_results/matrix/{sample}/iced/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_iced.matrix.biases",
        original="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/iced.matrix.biases.original",
        bias_vec="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/iced.matrix.biases",
        features_sorted="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/restriction_enzyme.fithichip.bed.sorted",
    log:
        "logs/fithichip_bias_ice.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        sed 's/nan/0/g' {params.ice} > {params.original} 2>> {log}
        paste {params.bin_intervals} {params.original} | cut -f1,2,3,5 > {params.bias_vec} 2>> {log}
        Rscript {script_dir}/../scripts/fithichip/BiasCalc.r --CoverageFile {input.coverage} --BiasFile {params.bias_vec} --OutFile {output.coverage_bias} 2>> {log}
        #Rscript {script_dir}/../scripts/fithichip/BiasCalc.r --CoverageFile {input.coverage} --BiasFile {params.bias_vec} --OutFile {output.coverage_bias}.temp 2>> {log}
        #sort -k1,1 -k2,2n {params.features} > {params.features_sorted}
        #awk 'NR>1' {output.coverage_bias}.temp | sort -k1,1 -k2,2n > {output.coverage_bias}.sorted 2>> {log}
        #bedtools map -c 4,5 -o mean -null '0' -a {output.coverage_bias}.sorted -b {params.features_sorted} > {output.coverage_bias}.temp1 2>> {log}
        #bedtools map -a {output.coverage_bias}.temp1 -b {params.features_sorted} -c 4 -o count -null '0' > {output.coverage_bias} 2>> {log}
        #rm {output.coverage_bias}.temp {output.coverage_bias}.sorted {output.coverage_bias}.temp1 2>> {log}
        """

rule fithichip_generate_features:
    input:
        coverage_bias = lambda wildcards: (
            f"fithichip/{wildcards.sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.coverage.bed"
            if fithichip_BiasType == 1
            else f"fithichip/{wildcards.sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.ice.bed"
        )
    output:
        features="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.allbins.genomefeature.bed",
    shell:
        """
        awk 'NR>1 {{print $0"\\t0\\t0\\t0"}}' {input.coverage_bias} | sort -k1,1 -k2,2n > {output.features}
        """

rule fithichip_interactions_with_features:
    input:
        interactions="fithichip/{sample}/Matrix_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.cis.interactions.DistThr.bed",
        features="fithichip/{sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.allbins.genomefeature.bed",
    output:
        interactions_features="fithichip/{sample}/Interactions_BinSize"+str(bin_size)+"/Interactions.{sample}.bed",
    params:
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
    log:
        "logs/fithichip_interactions_with_features.{sample}.log",
    conda:
        "envs/fithichip.yaml"  # Assumes an environment with R
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/Significance_Features.r -I {input.interactions} -E {input.features} -O {output.interactions_features} -C {params.chromosome_sizes} > {log} 2>&1
        """

rule fithichip_sort_interactions:
    input:
        interactions_features="fithichip/{sample}/Interactions_BinSize"+str(bin_size)+"/Interactions.{sample}.bed",
    output:
        sorted_interactions="fithichip/{sample}/Interactions_BinSize"+str(bin_size)+"/Interactions.{sample}.sortedGenDist.bed",
        peak_interactions="fithichip/{sample}/Interactions_BinSize"+str(bin_size)+"/Interactions.{sample}.sortedGenDist.peak.bed",
    params:
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
    log:
        "logs/fithichip_sort_interactions.{sample}.log",
    conda:"envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/Interaction_Sort_Genomic_Distance.r {input.interactions_features} {output.sorted_interactions} 0 {params.chromosome_sizes} 7 > {log} 2>&1
        awk '((NR==1) || ($9==1) || ($15==1))' {output.sorted_interactions} > {output.peak_interactions}  2>> {log}
        """

rule fithichip_significance_peak:
    input:
        sorted_interactions="fithichip/{sample}/Interactions_BinSize"+str(bin_size)+"/Interactions.{sample}.sortedGenDist.peak.bed",
        coverage_bias = lambda wildcards: (
            f"fithichip/{wildcards.sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.coverage.bed"
            if fithichip_BiasType == 1
            else f"fithichip/{wildcards.sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.ice.bed"
        )
    output:
        interactions="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.peaktoall.bedpe",
        sig_interactions="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.peaktoall.significant.bedpe",
    params:
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
        outdir="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size),
        qvalue=fdr,
        nbins=200,
        int_type=3,
        bias_corr=1,
        bias_type=fithichip_BiasType,
        bin_size=bin_size,
    log:
        "logs/fithichip_significance_peak.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/FitHiC_SigInt.r --InpFile {input.sorted_interactions} --headerInp \
            --OutFile {output.interactions} \
            --CoverageFile {input.coverage_bias} --BinSize {bin_size} --P2P 0 --IntType {params.int_type} \
            --BiasCorr {params.bias_corr} --BiasType {params.bias_type} --BiasLowThr 0.2 --BiasHighThr 5 \
            --UseNonzeroContacts 0 --Draw --cccol 7 --BiasFilt 0 --MultBias 0 --Resid 0 --EqOcc 1  > {log} 2>&1
        awk -F'\\t' -v q={params.qvalue} '{{if ((NR==1) || ($NF != "NA" && sprintf("%0.400f",$NF) < q)) {{print $0}}}}' {output.interactions} > {output.sig_interactions} 2>>{log}
        """

rule fithichip_significance_all:
    input:
        sorted_interactions="fithichip/{sample}/Interactions_BinSize"+str(bin_size)+"/Interactions.{sample}.sortedGenDist.bed",
        coverage_bias = lambda wildcards: (
            f"fithichip/{wildcards.sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.coverage.bed"
            if fithichip_BiasType == 1
            else f"fithichip/{wildcards.sample}/NormFeatures_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.bias.ice.bed"
        )
    output:
        interactions="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.alltoall.bedpe",
        sig_interactions="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.alltoall.significant.bedpe",
    params:
        chromosome_sizes=expand(genome_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
        outdir="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size),
        qvalue=fdr,
        nbins=200,
        int_type=4,
        bias_corr=1,
        bias_type=fithichip_BiasType,
        bin_size=bin_size,
    log:
        "logs/fithichip_significance_all.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
       Rscript {script_dir}/../scripts/fithichip/FitHiC_SigInt.r --InpFile {input.sorted_interactions} --headerInp \
            --OutFile {output.interactions} \
            --CoverageFile {input.coverage_bias} --BinSize {bin_size} --P2P 0 --IntType {params.int_type} \
            --BiasCorr {params.bias_corr} --BiasType {params.bias_type} --BiasLowThr 0.2 --BiasHighThr 5 \
            --UseNonzeroContacts 0 --Draw --cccol 7 --BiasFilt 0 --MultBias 0 --Resid 0 --EqOcc 1  > {log} 2>&1
        awk -F'\\t' -v q={params.qvalue} '{{if ((NR==1) || ($NF != "NA" && sprintf("%0.400f",$NF) < q)) {{print $0}}}}'  {output.interactions} > {output.sig_interactions} 2>>{log}
        """

rule fithichip_merge_interactions:
    input:
        sig_interactions_all="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.alltoall.significant.bedpe",
        sig_interactions_peak="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.peaktoall.significant.bedpe",
    output:
        merged_all="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.alltoall.significant.merged.bedpe",    
        merged_peak="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.peaktoall.significant.merged.bedpe",    
    params:
        bin_size=bin_size,
    log:
        "logs/fithichip_merge_interactions.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        python {script_dir}/../scripts/fithichip/CombineNearbyInteraction.py --InpFile {input.sig_interactions_peak} --OutFile {output.merged_peak} \
            --headerInp 1 --binsize {params.bin_size} --percent 100 --Neigh 2 > {log} 2>&1
        python {script_dir}/../scripts/fithichip/CombineNearbyInteraction.py --InpFile {input.sig_interactions_all} --OutFile {output.merged_all} \
            --headerInp 1 --binsize {params.bin_size} --percent 100 --Neigh 2 >> {log} 2>&1
        """

rule fithichip_summary_peak:
    input:
        sig_interactions_all="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.alltoall.significant.bedpe",
        sig_interactions_peak="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.peaktoall.significant.bedpe",
    output:
        distance_plot_all="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/Plots/FitHiChIP.{sample}.interactions.alltoall.significant.Dist_CC.png",
        distance_plot_peak="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/Plots/FitHiChIP.{sample}.interactions.peaktoall.significant.Dist_CC.png",
    log:
        "logs/fithichip_summary_peak.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/Distance_vs_CC.r --IntFile {input.sig_interactions_all} --OutFile {output.distance_plot_all} > {log} 2>&1 || true  
        Rscript {script_dir}/../scripts/fithichip/Distance_vs_CC.r --IntFile {input.sig_interactions_peak} --OutFile {output.distance_plot_peak} >> {log} 2>&1 || true        
       """

rule fithichip_summary_all:
    input:
        merged_all="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.alltoall.significant.merged.bedpe",    
        merged_peak="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/FitHiChIP.{sample}.interactions.peaktoall.significant.merged.bedpe",    
    output:
        merge_plot_all="fithichip/{sample}/FitHiChIP_alltoall_BinSize"+str(bin_size)+"/Plots/FitHiChIP.{sample}.interactions.alltoall.significant.merged.Dist_CC.png",
        merge_plot_peak="fithichip/{sample}/FitHiChIP_peaktoall_BinSize"+str(bin_size)+"/Plots/FitHiChIP.{sample}.interactions.peaktoall.significant.merged.Dist_CC.png",
    log:
        "logs/fithichip_summary_all.{sample}.log",
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/fithichip/Distance_vs_CC.r --IntFile {input.merged_all} --OutFile {output.merge_plot_all} > {log} 2>&1 || true
        Rscript {script_dir}/../scripts/fithichip/Distance_vs_CC.r --IntFile {input.merged_peak} --OutFile {output.merge_plot_peak} >> {log} 2>&1 || true
        """

#############################################
# --- MAPS Pipeline Rules ---
#############################################

rule feather_filter:
    input:
        fastqs=expand("processedFastqs/{{sample}}/{{sample}}_{read}.fastq.gz", read=reads),
    output:
        bam="maps_feather/{sample}/{sample}.paired.srtn.rmdup.bam",
    params:
        out_dir="maps_feather/{sample}",
        chromap_index=genome_dir + "/organisms/" + organism + "/ChromapIndex/" + organism + ".index",
        chromap_index_fa=genome_dir + "/organisms/" + organism + "/" + organism + ".fa"
    log:
        "logs/feather_filter.{sample}.log"
    benchmark:
        "benchmarks/feather_filter.{sample}.benchmark.txt"
    conda:
        "envs/maps.yaml"
    threads: threads
    shell:
        """
        mkdir -p {params.out_dir}
        python {script_dir}/../scripts/feather_filter_all.py \
            -o {params.out_dir} \
            -p {wildcards.sample} \
            -f1 {input.fastqs[0]} \
            -f2 {input.fastqs[1]} \
            -b {params.chromap_index} \
            -fa {params.chromap_index_fa} \
            -q {mapq} \
            -t {threads} \
            -d {optical_duplicate_distance} 2> {log}
        """

rule feather_split:
    input:
        bam="maps_feather/{sample}/{sample}.paired.srtn.rmdup.bam",
        peaks="peaks/{sample}_peaks.narrowPeak",
    output:
        bedpe="maps_feather/{sample}/{sample}.long.intra.bedpe",
        qc="maps_feather/{sample}/{sample}.feather.qc",
    params:
        out_dir="maps_feather/{sample}",
        peaks="peaks/{sample}_peaks.narrowPeak" if macs2_peaks == "NA" else macs2_peaks,
    log:
        "logs/feather_split.{sample}.log"
    benchmark:
        "benchmarks/feather_split.{sample}.benchmark.txt"
    conda:
        "envs/maps.yaml"
    threads: threads
    shell:
        """
        python {script_dir}/../scripts/feather_split_all.py \
            -o {params.out_dir} \
            -p {wildcards.sample} \
            -s {input.bam} \
            -l {length_cutoff} \
            -c True \
            -j {generate_hic} \
            -d {downsample_size} \
            -a {params.peaks} 2> {log}
        """

rule feather_process_qc:
    input:
        qc="maps_feather/{sample}/{sample}.feather.qc",
    output:
        qc="qc/{sample}_maps_specific_qc.tsv"
    params:
        log_file="maps_feather/{sample}/tempfiles/{sample}_R1.fastq.gz.sam.log"
    log:
        "logs/feather_process_qc.{sample}.log"
    benchmark:
        "benchmarks/feather_process_qc.{sample}.benchmark.txt"
    conda:
        "envs/maps.yaml"
    shell:
        """
            tail -n 12 {params.log_file} | sed -n '1p;2p;3p;4p;11p' | sed 's/\\.//g' | sed 's/: /\\t/g' | awk -F"\\t" 'BEGIN{{OFS="\\t"}} NR==1 {{total=$2/2; print $1, $2/2, "100.0%"}} NR>1 {{val=$2/2; pct=(val/total)*100; printf "%s\\t%d\\t%.2f%%\\n", $1, val, pct}}' > {output.qc} 2> {log}
            cat {input.qc} >> {output.qc} 2>> {log}
        """

# Rule to make MAPS runfile
rule maps_make_runfile:
    input:
        qc="maps_feather/{sample}/{sample}.long.intra.bedpe",
        fragment = expand(genome_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.organism.txt", 
                         organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size),
        peaks="peaks/{sample}_peaks.narrowPeak",
    output:
        "maps/{sample}/maps_{sample}.maps"
    params:
        out_dir = "maps/{sample}",
        feather_dir = "maps_feather/{sample}",
        chr_count = expand(genome_dir+"/organisms/{organism}/{organism}.autosome.count.txt", organism=organism),  # Dependency on calculated count
        peaks="peaks/{sample}_peaks.narrowPeak" if macs2_peaks == "NA" else macs2_peaks,
    log:
        "logs/maps_make_runfile.{sample}.log"
    benchmark:
        "benchmarks/maps_make_runfile.{sample}.benchmark.txt"
    conda:
        "envs/maps.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        python {script_dir}/../scripts/make_maps_runfile.py \
            --BINNING_RANGE {binning_range} \
            {wildcards.sample} \
            {params.out_dir}/ \
            {params.peaks} \
            {input.fragment} \
            {params.feather_dir}/ \
            {params.feather_dir}/ \
            {bin_size} \
            {params.chr_count} \
            {params.out_dir}/ \
            {maps_sex_chroms} \
            > {log} 2>&1
        """

rule maps_run_maps:
    input:
        "maps/{sample}/maps_{sample}.maps"
    output:
        "maps/{sample}/{sample}.maps.qc"
    log:
        "logs/maps_run_maps.{sample}.log"
    benchmark:
        "benchmarks/maps_run_maps.{sample}.benchmark.txt"
    conda:
        "envs/maps.yaml"
    shell:
        """
        python {script_dir}/../scripts/MAPS.py {input} > {log} 2>&1
        """

rule maps_regression:
    input:
        maps="maps/{sample}/{sample}.maps.qc",
        chr_count = expand(genome_dir+"/organisms/{organism}/{organism}.autosome.count.txt", organism=organism)  # Dependency on calculated count
    output:
        "maps/{sample}/summary.{sample}."+str(resolution)+"k.txt"
    params:
        out_dir="maps/{sample}",
    log:
        "logs/maps_regression.{sample}.log"
    benchmark:
        "benchmarks/maps_regression.{sample}.benchmark.txt"
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/MAPS_regression_and_peak_caller.r \
            {params.out_dir}/ \
            {wildcards.sample}.{resolution}k \
            {bin_size} \
            {input.chr_count} \
            {maps_sex_chroms} \
            {maps_model} \
            {maps_count_cutoff} \
            {maps_ratio_cutoff} \
            {fdr} \
            None > {log} 2>&1
        """

rule maps_peak_formatting:
    input:
        qc_report="qc/{sample}_maps_specific_qc.tsv",
        summary="maps/{sample}/summary.{sample}."+str(resolution)+"k.txt"
    output:
        "maps/{sample}/Maps.{sample}."+str(resolution)+"k.sig3Dinteractions.bedpe"
    params:
        out_dir="maps/{sample}",
    log:
        "logs/maps_peak_formatting.{sample}.log"
    benchmark:
        "benchmarks/maps_peak_formatting.{sample}.benchmark.txt"
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/MAPS_peak_formatting.r \
            {params.out_dir}/ \
            {wildcards.sample}.{resolution}k \
            {fdr} \
            {bin_size} > {log} 2>&1
        """

###############################################################################
# --- Format matrix ---
###############################################################################
# Rule to sort and index BAM files
rule sort_index_bam:
    input:
        bam="hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.bam",
    output:
        sorted_bam=temp("hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.sorted.bam"),
        index=temp("hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.sorted.bam.bai")
    conda:
        "envs/hicpro.yaml"
    threads: threads
    shell:
        """
        mkdir -p $(dirname {output.sorted_bam})
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam} &&
        samtools index -@ {threads} {output.sorted_bam}
        """

# Rule to generate bigWig files
rule generate_bigwig:
    input:
        sorted_bam="hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.sorted.bam"
    output:
        bigwig="matrices/{sample}.bwt2pairs.bw"
    params:
        out_dir = "hicpro/{sample}/bowtie_results/bwt2/{sample}",
        genome_size = lambda wildcards: sum(
            int(line.strip().split()[1])
            for line in open(f"{genome_dir}/organisms/{organism}/{organism}.chrom.sizes")
            if line.strip()
        ),
    conda:
        "envs/deeptools.yaml"
    threads: threads
    shell:
        """
        bamCoverage -b {input.sorted_bam} -o {output.bigwig} \
            --binSize {BINSIZE} \
            --effectiveGenomeSize {params.genome_size} \
            -p {threads}
        """

# Rule to convert HiC-Pro matrices to cool format
rule convert_to_cool:
    input:
        matrix="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+".matrix",
        bed = "hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_abs.bed",
    output:
        cool="matrices/{sample}_"+str(bin_size)+".cool"
    conda:
        "envs/hicexplorer.yaml"
    threads: threads
    shell:
        """
        hicConvertFormat --matrices {input.matrix} \
            --outFileName {output.cool} \
            --inputFormat hicpro \
            --outputFormat cool \
            --bedFileHicpro {input.bed}
        """

# Rule to convert HiC-Pro matrices to h5 format
rule convert_to_h5:
    input:
        matrix="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+".matrix",
        bed="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_abs.bed"
    output:
        h5="matrices/{sample}_"+str(bin_size)+".h5"
    conda:
        "envs/hicexplorer.yaml"
    threads: threads
    shell:
        """
        hicConvertFormat --matrices {input.matrix} \
            --outFileName {output.h5} \
            --inputFormat hicpro \
            --outputFormat h5 \
            --bedFileHicpro {input.bed}
        """

###############################################################################
# --- QC results quantile---
###############################################################################
# Rule to process quantiles
rule process_quantiles:
    input:
        matrix="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+".matrix",
        bed="hicpro/{sample}/hic_results/matrix/{sample}/raw/"+str(bin_size)+"/{sample}_"+str(bin_size)+"_abs.bed",
        cool="matrices/{sample}_"+str(bin_size)+".cool",
        h5="matrices/{sample}_"+str(bin_size)+".h5"
    output:
        nonzero="quantiles/{sample}/hichip_matrix/quantiles_nonzero.csv",
        all_csv="quantiles/{sample}/hichip_matrix/quantiles_all.csv"
    params:
        out_dir="quantiles/{sample}/hichip_matrix",
        matrix="matrices/{sample}_"+str(bin_size)+".matrix",
        bed="matrices/{sample}_"+str(bin_size)+"_abs.bed",
        valid_pairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
        valid_pairs_link="matrices/{sample}.allValidPairs",
    conda:
        "envs/cooler.yaml"
    shell:
        """
        ln -sf ../{params.valid_pairs} {params.valid_pairs_link}
        ln -sf ../{input.matrix} {params.matrix}
        ln -sf ../{input.bed} {params.bed}
        python {script_dir}/../scripts/qc_hichip_quantiles.py --input {input.cool} --output_dir {params.out_dir}
        """

# Rule to generate QC plots
rule generate_qc_plots1:
    input:
        nonzero_files=expand("quantiles/{sample}/hichip_matrix/quantiles_nonzero.csv", sample=samples),
        all_files=expand("quantiles/{sample}/hichip_matrix/quantiles_all.csv", sample=samples)
    output:
        plot="quantiles/quantile_dot_plot_hichip.pdf"
    params:
        outdir="quantiles",
        files="quantiles/*/hichip_matrix/quantiles_nonzero.csv",
        name="hichip",
        all_files="quantiles/*/hichip_matrix/quantiles_all.csv"
    log:
        os.path.join("logs", "quantiles1.log")
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/qc_plot_quantiles.R {params.outdir} "{params.files}" {params.name} > {log} 2>&1
        """

rule generate_qc_plots2:
    input:
        nonzero_files=expand("quantiles/{sample}/hichip_matrix/quantiles_nonzero.csv", sample=samples),
        all_files=expand("quantiles/{sample}/hichip_matrix/quantiles_all.csv", sample=samples)
    output:
        plot="quantiles/quantile_dot_plot_hichip_over1.pdf"
    params:
        outdir="quantiles",
        files="quantiles/*/hichip_matrix/quantiles_all.csv",
        name="hichip_over1",
    log:
        os.path.join("logs", "quantiles2.log")
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/qc_plot_quantiles.R {params.outdir} "{params.files}" {params.name} > {log} 2>&1
        """

###############################################################################
# --- QC results multiQC reads---
###############################################################################
rule fastqc:
    input:
        fastq="originalFastqs/{sample}_{read}.fastq.gz" #if int(downsample_size) == 0 else "downsampledFastqs/{sample}_{read}.fastq.gz"
    output:
        "fastqc/{sample}_{read}_fastqc.html"
    log:
        "logs/fastqc.{sample}_{read}.log"
    benchmark:
        "benchmarks/fastqc.{sample}_{read}.txt"
    threads: threads
    conda:
        "envs/basic.yaml"
    shell:
        """
        fastqc -t {threads} {input.fastq} -o fastqc 2>> {log}
        """

def get_optional_inputs(wildcards):
    """Return a list of optional input files based on pipeline configuration."""
    possible_inputs = {
        "hicdcplus_loops": f"hicdcplus/{wildcards.sample}/Hicdcplus.{wildcards.sample}.significant_interactions.bedpe",
        "hichipper_loops": f"hichipper/Hichipper.{wildcards.sample}.significant.interactions.bedpe",
        "fithichip_loops": f"fithichip/{wildcards.sample}/FitHiChIP_peaktoall_BinSize{bin_size}/FitHiChIP.{wildcards.sample}.interactions.peaktoall.significant.bedpe",
        "maps_loops": f"maps/{wildcards.sample}/Maps.{wildcards.sample}.{resolution}k.sig3Dinteractions.bedpe",
    }
    pipeline = config.get("pipeline", "All")
    active_inputs = []
    if pipeline == "All":
        # Return all possible inputs when pipeline is "All"
        active_inputs = list(possible_inputs.values())
    else:
        # Return the input corresponding to the pipeline (case-insensitive match)
        for key, path in possible_inputs.items():
            if pipeline.upper() == key.split("_")[0].upper():
                active_inputs.append(path)
                break  # Only include the matching pipeline’s input
    return active_inputs

rule hichip_qc:
    input:
        valid_pairs="hicpro/{sample}/hic_results/data/{sample}/{sample}.allValidPairs",
        peaks="peaks/{sample}_peaks.narrowPeak",
        optional_inputs=get_optional_inputs,  # Dynamically determine optional files
    output:
        qc="qc/{sample}_hichip_qc.tsv"
    params:
        mapping="hicpro/{sample}/bowtie_results/bwt2/{sample}/{sample}_"+organism+".bwt2pairs.pairstat",
        interactions="hicpro/{sample}/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat",
    log:
        "logs/hichip_qc.{sample}.log"
    conda:
        "envs/r.yaml"
    shell:
        """
        mkdir -p significant_interactions && cp {input.optional_inputs} significant_interactions/
        Rscript {script_dir}/../scripts/qc_all.R \
            {wildcards.sample} \
            {output.qc} \
            {params.mapping} \
            {params.interactions} \
            {input.peaks} \
            {input.optional_inputs} > {log} 2>&1
        """
                
rule multiqc:
    input:
        fastqc=expand("fastqc/{sample}_{read}_fastqc.html", sample=samples, read=reads),
        hichip_qc=expand("qc/{sample}_hichip_qc.tsv", sample=samples),
        maps_qc=expand("qc/{sample}_maps_specific_qc.tsv", sample=samples) if pipeline in ["All", "Maps"] else []
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/r.yaml"
    shell:
        """
        multiqc . --ignore "*flagstat*" --force --outdir qc/ --filename multiqc_report.html > {log} 2>&1
        """

# Rule to generate valid pairs barplot
rule generate_valid_pairs_barplot:
    input:
        hichip_qc=expand("qc/{sample}_hichip_qc.tsv", sample=samples),
    output:
        barplot="qc/summary_hichip_valid_pairs_barplot.pdf"
    params:
        qc_files="qc/*hichip_qc.tsv",
        qc_dir="qc",
        pattern="hichip_qc.tsv",
        out_dir="qc",
        out_prefix="hichip_valid_pairs",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/qc_hichip_barplot.R {params.qc_dir} {params.pattern} Valid_Pairs {params.out_dir} {params.out_prefix}
        """

# Rule to generate FRiP barplot
rule generate_frip_barplot:
    input:
        hichip_qc=expand("qc/{sample}_hichip_qc.tsv", sample=samples),
    output:
        barplot="qc/summary_hichip_frip_barplot.pdf"
    params:
        qc_files="qc/*maps_specific_qc.tsv",
        qc_dir="qc",
        pattern="maps_specific_qc.tsv",
        out_dir="qc",
        out_prefix="hichip_frip",
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/qc_hichip_barplot.R {params.qc_dir} {params.pattern} FRiP {params.out_dir} {params.out_prefix}
        """



#############################################
# --- Consensus interaction calling ---
# Keep interactions supported by >= min_support callers (default: 2 of 4).
# Outputs:
#   consensus_interactions/{sample}.consensus_interactions.bedpe        <- final list
#   consensus_interactions/{sample}.consensus_interactions.all_callers.bedpe  <- all + support counts
#############################################

rule consensus_interactions:
    input:
        # Each caller file is optional: if a caller was not run, pass an empty
        # string and the Python script will skip it gracefully.
        fithichip  = lambda wc: f"significant_interactions/FitHiChIP.{wc.sample}.interactions.peaktoall.significant.bedpe"
                                 if pipeline in ["All", "Fithichip"] else [],
        hicdcplus  = lambda wc: f"significant_interactions/Hicdcplus.{wc.sample}.significant_interactions.bedpe"
                                 if pipeline in ["All", "Hicdcplus"] else [],
        hichipper  = lambda wc: f"significant_interactions/Hichipper.{wc.sample}.significant.interactions.bedpe"
                                 if pipeline in ["All", "Hichipper"] else [],
        maps       = lambda wc: f"significant_interactions/Maps.{wc.sample}.{resolution}k.sig3Dinteractions.bedpe"
                                 if pipeline in ["All", "Maps"] else [],
    output:
        consensus  = "consensus_interactions/{sample}.consensus_interactions.bedpe",
    params:
        binsize    = bin_size,
        min_support= 2 # #config.get("consensus_min_support", 2),  # set in config.yaml
    log:
        "logs/consensus_interactions.{sample}.log"
    benchmark:
        "benchmarks/consensus_interactions.{sample}.txt"
    conda:
        "envs/basic.yaml"   # needs python + pandas (already in most envs)
    shell:
        """
        mkdir -p consensus_interactions

        python {script_dir}/../scripts/consensus_interactions.py \
            {("--fithichip "  + input.fithichip  ) if input.fithichip  else ""} \
            {("--hicdcplus "  + input.hicdcplus  ) if input.hicdcplus  else ""} \
            {("--hichipper "  + input.hichipper  ) if input.hichipper  else ""} \
            {("--maps "       + input.maps       ) if input.maps        else ""} \
            --binsize     {params.binsize}      \
            --min-support {params.min_support}  \
            --out         {output.consensus}    \
            > {log} 2>&1
        """














#############################################
# --- Differential analysis using hicdcplus Pipeline Rules ---
#############################################
# Rule to identify differential interactions (conditional on samples_comparison)
rule hicdc_differential_interactions:
    input:
        control_gz_files=glob.glob("hicdcplus/{samples_comparison[0]}*/Hicdcplus.{samples_comparison[0]}*significant_interactions.bedpe.gz"),
        test_gz_files=glob.glob("hicdcplus/{samples_comparison[1]}*/Hicdcplus.{samples_comparison[1]}*significant_interactions.bedpe.gz")
    output:
        bedpe="differential_interactions_{samples_comparison[0]}_vs_{samples_comparison[1]}_hicdcplus/dispersionplot.pdf",
    params:
        outdir="differential_interactions_{samples_comparison[0]}_vs_{samples_comparison[1]}_hicdcplus/",
        fdr=fdr,  
        control=samples_comparison[0],
        test=samples_comparison[1],
    log:
        "logs/hicdc_differential_interactions.{samples_comparison[0]}_vs_{samples_comparison[1]}.log"
    benchmark:
        "benchmarks/hicdc_differential_interactions.{samples_comparison[0]}_vs_{samples_comparison[1]}.txt"
    threads: 1
    conda:
        "envs/hicdcplus.yaml"
    shell:
        """
        Rscript {script_dir}/../scripts/hicdc_differential_interactions.R \
          --control {params.control} \
          --test {params.test} \
          --control_gz_files "{input.control_gz_files_str}" \
          --test_gz_files "{input.test_gz_files_str}" \
          --output_dir {params.outdir} \
          --fdr {params.fdr} > {log} 2>&1
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



