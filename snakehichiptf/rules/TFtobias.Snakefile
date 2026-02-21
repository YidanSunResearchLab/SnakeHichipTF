import os
import yaml
from snakemake.utils import validate
shell.executable("bash") ## Force bash as the shell for all rules

#############################################
# Configuration
#############################################
script_dir = os.path.dirname(workflow.snakefile)

# Load configuration
organism = config["organism"]
genome_dir = config["genome_dir"]
atacseq_bam = config["atacseq_bam"].split()  # Split space-separated BAM files
enhancer_bed = config["enhancer_bed"].split()  # Split space-separated BEDPE files
samples = config["samples"]
motifs = config["motifs"]
atacorrect = config["atacorrect"]
footprinting = config["footprinting"]
plot_aggregate = config["plot_aggregate"]
plot_changes = config["plot_changes"]
bindetect = config["bindetect"]
threads = config["threads"]  # Added threads from config
output_dir = os.path.join("tobias", "_".join(samples))
default_motifs = os.path.join(script_dir, "../organisms", "all_motifs.txt")

# Default rule
rule all:
    input:
        # BINDetect output (single directory)
        os.path.join(output_dir, "bindetect","bindetect_results.txt"),
        # Join bound outputs
        expand(os.path.join(output_dir, "all_{sample}_bound.bed"), sample=samples),
        # PlotAggregate outputs
        expand(os.path.join(output_dir, "plots", "{sample}_aggregate.pdf"), sample=samples),
        # PlotChanges output
        #os.path.join(output_dir, "plots", "changes_{control}_vs_{test}.pdf").format(control=control, test=test)

# 1. Link BAM files
rule link_bam:
    input:
        bam = lambda wildcards: atacseq_bam[samples.index(wildcards.sample)]
    output:
        bam = os.path.join(output_dir, "{sample}.bam")
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "link_bam_{sample}.log")
    benchmark: os.path.join(output_dir, "benchmarks", "link_bam_{sample}.txt")
    shell:
        "ln -sf {input.bam} {output.bam} 2> {log}"

# 2. Merge BEDPE files and convert to BED (simplified)
rule merge_enhancer_beds:
    input:
        bedpe_files = enhancer_bed
    output:
        merged = os.path.join(output_dir, "merged_enhancers.bed")
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "merge_enhancer_beds.log")
    benchmark: os.path.join(output_dir, "benchmarks", "merge_enhancer_beds.txt")
    shell:
        r"""
        cat {input.bedpe_files} | \
        awk -F'\t' 'NF>=6 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $5~/^[0-9]+$/ && $6~/^[0-9]+$/ \
            {{print $1"\t"$2"\t"$3"\t"$1"_"$2"\t0\t+" "\n" $4"\t"$5"\t"$6"\t"$4"_"$5"\t0\t+"}}' | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i - | \
        awk -F'\t' '{{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"\t0\t+"}}' \
        > {output.merged} 2> {log}
        """
      
# 3. Format motifs
rule format_motifs:
    output:
        joined = os.path.join("tobias", "motifs", "all_motifs.txt")
    params:
        motifdir = os.path.join("tobias", "motifs", "individual"),
        default_motifs = default_motifs,
        motifs = motifs,
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "format_motifs.log")
    benchmark: os.path.join(output_dir, "benchmarks", "format_motifs.txt")
    localrule: True
    shell:
        """
        mkdir -p {params.motifdir} && \
        if [ "{params.motifs}" = "NA" ]; then
            cp {params.default_motifs} {output.joined} 2> {log}
        else
            echo "Copying from {params.motifs} to {params.motifdir}" &>> {log} && \
            cp {params.motifs}/* {params.motifdir} 2>> {log} && \
            format=$(for f in {params.motifdir}/*; do grep -oE '\.(jaspar|pfm|meme)$' <<< "$f" | sed 's/\.//'; break; done) && echo "$format" &>> {log}
            if [ -z "$format" ]; then
                echo "Error: Unknown motif format in {params.motifdir}" &>> {log}
                exit 1
            else
                TOBIAS FormatMotifs --input {params.motifdir} --format "$format" --task join --output {output.joined} &>> {log}
            fi
        fi
        
        """
        
# 4. TOBIAS ATACorrect
rule tobias_atacorrect:
    input:
        bam = os.path.join(output_dir, "{sample}.bam"),
        genome = expand(genome_dir + "/organisms/{organism}/{organism}.fa", organism=organism),
        peaks = rules.merge_enhancer_beds.output.merged
    output:
        corrected_bw = os.path.join(output_dir, "bias_correction", "{sample}_corrected.bw")
    params:
        extra = "" if atacorrect == "NA" else atacorrect,
        outdir = os.path.join(output_dir, "bias_correction"),
        prefix = lambda wildcards: wildcards.sample,
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "atacorrect_{sample}.log")
    benchmark: os.path.join(output_dir, "benchmarks", "atacorrect_{sample}.txt")
    shell:
        "TOBIAS ATACorrect -b {input.bam} -g {input.genome} -p {input.peaks} {params.extra} "
        "--outdir {params.outdir} --prefix {params.prefix} --cores {threads} &>> {log}"

# 5. TOBIAS FootprintScores
rule tobias_footprint:
    input:
        corrected_bw = rules.tobias_atacorrect.output.corrected_bw,
        peaks = rules.merge_enhancer_beds.output.merged
    output:
        footprints_bw = os.path.join(output_dir, "footprinting", "{sample}_footprints.bw")
    params:
        extra = "" if footprinting == "NA" else footprinting
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "footprint_{sample}.log")
    benchmark: os.path.join(output_dir, "benchmarks", "footprint_{sample}.txt")
    shell:
        "TOBIAS FootprintScores -s {input.corrected_bw} -r {input.peaks} "
        "{params.extra} -o {output.footprints_bw} --cores {threads} &>> {log}"

# 6. TOBIAS BINDetect (across all samples)
def get_bindetect_params(samples, extra):
    extra = "" if extra == "NA" else extra
    if len(samples) < 3:
        return f"--cond_names {' '.join(samples)} {extra}"
    else:
        return f"--time-series {extra}"

rule tobias_bindetect:
    input:
        footprints = expand(os.path.join(output_dir, "footprinting", "{sample}_footprints.bw"), sample=samples),
        motifs = os.path.join("tobias", "motifs", "all_motifs.txt"),
        genome = expand(genome_dir + "/organisms/{organism}/{organism}.fa", organism=organism),
        peaks = rules.merge_enhancer_beds.output.merged,
    output:
        result = os.path.join(output_dir, "bindetect", "bindetect_results.txt")
    params:
        outdir = os.path.join(output_dir, "bindetect"),
        bindetect_params = lambda wildcards, input: get_bindetect_params(samples, bindetect),
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "bindetect.log")
    benchmark: os.path.join(output_dir, "benchmarks", "bindetect.txt")
    shell:
        "TOBIAS BINDetect --signals {input.footprints} --motifs {input.motifs} --genome {input.genome} "
        "--peaks {input.peaks} --outdir {params.outdir} --cores {threads} {params.bindetect_params} &>> {log}"

# 7. Checkpoint to get TF IDs
checkpoint get_TF_ids:
    input:
        bindetect = rules.tobias_bindetect.output.result
    output:
        tf_list = os.path.join(output_dir, "bindetect", "tfs_{sample}.txt")
    params:
        bindetect = os.path.join(output_dir, "bindetect")
    threads: 1
    log: os.path.join(output_dir, "logs", "get_tf_ids_{sample}.log")
    benchmark: os.path.join(output_dir, "benchmarks", "get_tf_ids_{sample}.txt")
    run:
        bindetect_dir = os.path.dirname(os.path.dirname(params.bindetect))  # Get bindetect directory
        # Corrected glob_wildcards without sample keyword
        wildcard_results = glob_wildcards(os.path.join(bindetect_dir, "{TF}", "beds", "{TF}_{sample}_bound.bed"))
        # Filter TFs for the current sample
        tf_ids = [tf for tf, sample in zip(wildcard_results.TF, wildcard_results.sample) if sample == wildcards.sample]
        with open(output.tf_list, "w") as f:
            f.write("\n".join(tf_ids))

def get_TF_ids(wildcards):
    tf_file = checkpoints.get_TF_ids.get(**wildcards).output.tf_list
    with open(tf_file, "r") as f:
        return [line.strip() for line in f if line.strip()]

# 8. Join bound estimates per sample
rule join_bound:
    input:
        tf_list = os.path.join(output_dir, "bindetect", "tfs_{sample}.txt"),
        beds = lambda wildcards: expand(os.path.join(output_dir, "bindetect", "{TF}", "beds", "{TF}_{sample}_bound.bed"), 
                                        TF=get_TF_ids(wildcards), sample=wildcards.sample)
    output:
        unsorted = temp(os.path.join(output_dir, "all_{sample}_bound.tmp")),
        final = os.path.join(output_dir, "all_{sample}_bound.bed")
    threads: threads
    log: os.path.join(output_dir, "logs", "join_bound_{sample}.log")
    benchmark: os.path.join(output_dir, "benchmarks", "join_bound_{sample}.txt")
    run:
        n = 100  # chunks of 100 files
        shell("> {output.unsorted}")
        file_chunks = [input.beds[i:i+n] for i in range(0, len(input.beds), n)]
        for chunk in file_chunks:
            shell("cat {chunk} >> {output.unsorted}")
        shell("sort -k1,1 -k2,2n {output.unsorted} > {output.final}")
        #shell("igvtools index {output.final} 2>> {log}")
        
# 9. TOBIAS PlotAggregate
rule tobias_plot_aggregate:
    input:
        footprints = rules.tobias_footprint.output.footprints_bw,
        bound = rules.join_bound.output.final  # Use joined bound BED instead of bindetect directory
    output:
        plot = os.path.join(output_dir, "plots", "{sample}_aggregate.pdf")
    params:
        extra = "" if plot_aggregate == "NA" else plot_aggregate
    threads: threads
    conda: "envs/tobias.yaml"
    log: os.path.join(output_dir, "logs", "plot_aggregate_{sample}.log")
    benchmark: os.path.join(output_dir, "benchmarks", "plot_aggregate_{sample}.txt")
    shell:
        "TOBIAS PlotAggregate --signals {input.footprints} --TFBS {input.bound} "
        "--output {output.plot} --plot_boundaries --share_y both {params.extra} &>> {log}"

