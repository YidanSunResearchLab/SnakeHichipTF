import json
import os
import yaml
from snakemake.utils import validate
shell.executable("bash")  # Force bash as the shell for all rules

#############################################
# Configuration
#############################################
script_dir = os.path.dirname(workflow.snakefile)

# Load configuration from wrapper
organism = config["organism"]
genome_dir = config["genome_dir"]
enhancer_bed = config["enhancer_bedpe"]
input_type = config["input_type"]  # "bam" or "fragments_tsv"

# wrapper now passes native types (mode is a plain string; samples is json.dumps(list) in your wrapper)
# sample_map is still json.dumps(dict) in your wrapper, so keep json.loads for sample_map
sample_map = json.loads(config["sample_map"]) if isinstance(config["sample_map"], str) else config["sample_map"]
samples = json.loads(config["samples"]) if isinstance(config["samples"], str) else config["samples"]

# mode is already a plain string after your wrapper change
mode = config.get("mode", "bulk")  # "bulk" or "singlecell"

threads = int(config.get("threads", 30))

# IMPORTANT: main_dir == output_dir == snakemake working directory
# Do NOT introduce a separate main_dir; just use os.getcwd() everywhere.
output_dir = os.path.join("scprinter", "_".join(samples))
main_dir = output_dir
entity = "test"
project = "scPrinter_ATAC"


# ------------------------------------------------------------
# Build a samples_cfg dict compatible with your old logic
# samples_cfg[s]["bam"] / ["fragments"]
# ------------------------------------------------------------
samples_cfg = {}
for s in samples:
    inp = sample_map[s]["input"]
    if input_type == "bam":
        samples_cfg[s] = {"bam": inp, "fragments": ""}
    else:
        samples_cfg[s] = {"bam": "", "fragments": inp}

rule all:
    input:
        f"{main_dir}/seq2print/scprinter_done.txt",
        f"{main_dir}/seq2print/TFBS_scores.csv",
        f"{main_dir}/seq2print/TFBS_scatter.png",

rule merge_enhancer_beds:
    input:
        bedpe_files = enhancer_bed
    output:
        merged = os.path.join(output_dir, "merged_enhancers.bed")
    threads: 1
    conda: "envs/maps.yaml"
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
      

if input_type == "bam":
    rule bam_to_fragments:
        input:
            bam=lambda wc: samples_cfg[wc.sample]["bam"]
        output:
            fragments = f"{main_dir}/fragments/{{sample}}_fragments.tsv.gz",
            index     = f"{main_dir}/fragments/{{sample}}_fragments.tsv.gz.tbi",
        params:
            name_sorted = f"{main_dir}/fragments/{{sample}}.name_sorted.bam",
            bed         = f"{main_dir}/fragments/{{sample}}_fragments.bed",
            tsv         = f"{main_dir}/fragments/{{sample}}_fragments.tsv",
        threads: threads
        conda: "envs/maps.yaml"
        log: f"{main_dir}/logs/bam_to_fragments/{{sample}}.log"
        benchmark: f"{main_dir}/benchmarks/bam_to_fragments/{{sample}}.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p {main_dir}/fragments {main_dir}/logs {main_dir}/benchmarks/bam_to_fragments

            samtools sort -n -@ {threads} -o {params.name_sorted} {input.bam}
            bedtools bamtobed -i {params.name_sorted} -bedpe > {params.bed}

            awk -v s="{wildcards.sample}" '$1 == $4 && $2 < $6 {{print $1"\t"$2"\t"$6"\t"s}}' {params.bed} \
                | sort -k1,1 -k2,2n > {params.tsv}

            bgzip -f {params.tsv}
            tabix -f -p bed {output.fragments} 2>> {log}
            """


if input_type == "fragments_tsv":
    rule link_fragments:
        input:
            src=lambda wc: samples_cfg[wc.sample]["fragments"]
        output:
            dest = f"{main_dir}/fragments/{{sample}}_fragments.tsv.gz"
        threads: 1
        log: f"{main_dir}/logs/link_fragments/{{sample}}.log"
        benchmark: f"{main_dir}/benchmarks/link_fragments/{{sample}}.txt"
        shell:
            r"""
            set -euo pipefail
            mkdir -p {main_dir}/fragments {main_dir}/logs/link_fragments {main_dir}/benchmarks/link_fragments
            ln -sf {input.src} {output.dest} 2> {log}
            """


if mode == "bulk":
    rule run_scprinter:
        input:
            sample = expand(f"{main_dir}/fragments/{{sample}}_fragments.tsv.gz", sample=samples),
            merged = os.path.join(output_dir, "merged_enhancers.bed")
        output:
            flag        = f"{main_dir}/seq2print/scprinter_done.txt",
            tfbs_scores = f"{main_dir}/seq2print/TFBS_scores.csv",
            tfbs_scatter= f"{main_dir}/seq2print/TFBS_scatter.png",
        threads: threads
        log:
            out = f"{main_dir}/logs/run_scprinter.out",
            err = f"{main_dir}/logs/run_scprinter.err"
        params:
            main_dir       = os.path.abspath(main_dir),
            entity         = entity,
            project        = project,
            footprints_arg = f"--footprints_bed {config['footprints_bed']}" if config.get("footprints_bed") else "",
            regions_arg    = f"--regions_bed {os.path.join(output_dir, 'merged_enhancers.bed')}" if config.get("enhancer_bed") else "",
            n_regions      = 10,
        conda:
            "envs/scprinter_gpu.yaml"
        benchmark: os.path.join(output_dir, "benchmarks", "run_scprinter.txt")
        shell:
            r"""
            set -euo pipefail

            out_file="{log.out}"
            err_file="{log.err}"

            i=1
            while [ -e "$out_file" ] || [ -e "$err_file" ]; do
                out_file="{log.out}_$i"
                err_file="{log.err}_$i"
                i=$((i+1))
            done

            echo "[INFO] Using log files:"
            echo "  stdout -> $out_file"
            echo "  stderr -> $err_file"

            exec > >(tee "$out_file") 2> >(tee "$err_file" >&2)
            echo "===== $(date) ====="

            export LD_LIBRARY_PATH="${{CONDA_PREFIX}}/lib:$LD_LIBRARY_PATH"
            export NUMBA_DISABLE_JIT=1
            echo "[DEBUG] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

            echo "[DEBUG] CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
            nvidia-smi || true

            python {script_dir}/../scripts/scprinter_main.py \
              --main_dir {params.main_dir} \
              --entity {params.entity} \
              --project {params.project} \
              {params.regions_arg} \
              {params.footprints_arg} \
              --n_regions {params.n_regions} \
              --resume

            touch {output.flag}
            """
else:
    rule run_scprinter:
        input:
            sample = expand(f"{main_dir}/fragments/{{sample}}_fragments.tsv.gz", sample=samples),
            merged = os.path.join(output_dir, "merged_enhancers.bed")
        output:
            flag         = f"{main_dir}/seq2print/scprinter_done.txt",
            tfbs_scores  = f"{main_dir}/seq2print/TFBS_scores.csv",
            tfbs_scatter = f"{main_dir}/seq2print/TFBS_scatter.png",
        log:
            out = f"{main_dir}/logs/run_scprinter.out",
            err = f"{main_dir}/logs/run_scprinter.err"
        threads: threads
        params:
            main_dir      = lambda wc: main_dir,
            device        = lambda wc: config.get("device", "cuda"),
            min_fragments = lambda wc: config.get("min_fragments", 1000),
            min_tsse      = lambda wc: config.get("min_tsse", 7),
            lora_model    = lambda wc: config.get("lora_model_path", ""),
            lora_config   = lambda wc: config.get("lora_config_path", ""),
            marker_dict   = lambda wc: json.dumps(config.get("marker_dict", {})),
        conda:
            "envs/scprinter_gpu.yaml"
        benchmark: os.path.join(output_dir, "benchmarks", "run_scprinter.txt")
        shell:
            r"""
            set -euo pipefail

            echo "[DEBUG] CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
            nvidia-smi || true

            export LD_LIBRARY_PATH="${{CONDA_PREFIX}}/lib:$LD_LIBRARY_PATH"
            echo "[DEBUG] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
            export NUMBA_DISABLE_JIT=1

            python {script_dir}/../scripts//singlecell_scprinter_main.py \
              --main_dir {params.main_dir} \
              --fragment_file {input[0]} \
              --device {params.device} \
              --min_fragments {params.min_fragments} \
              --min_tsse {params.min_tsse} \
              --resume \
              >> {log.out} 2>> {log.err}

            touch {output.flag}
            """





