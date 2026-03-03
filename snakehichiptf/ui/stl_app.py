import os
import shlex
import subprocess
import tempfile
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="SnakeHichipTF UI", layout="wide")
st.title("SnakeHichipTF — UI")
st.caption("Run Genomesetup / Snakehichip / Diffhichip / TFscprinter / TFtobias without the command line.")

# ---------- Execution helpers ----------

def render_cmd(cmd_list):
    st.code(" ".join(shlex.quote(c) for c in cmd_list))

def run_and_stream_local(cmd_list, workdir=None, env=None):
    st.write("**Command to run (local):**")
    render_cmd(cmd_list)

    log_box = st.empty()
    lines = []

    with subprocess.Popen(
        cmd_list,
        cwd=str(workdir) if workdir else None,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    ) as p:
        assert p.stdout is not None
        for line in p.stdout:
            lines.append(line.rstrip("\n"))
            log_box.code("\n".join(lines[-300:]))
        ret = p.wait()

    if ret == 0:
        st.success("Done.")
    else:
        st.error(f"Failed (exit code {ret}).")

def submit_slurm(cmd_list, job_name="snakehichip_ui", cpus=8, mem_gb=32, time="24:00:00", partition="", workdir=None):
    """
    Submit the command via sbatch. Good for long jobs on HPC.
    """
    st.write("**Command to submit (slurm):**")
    render_cmd(cmd_list)

    workdir = Path(workdir or os.getcwd()).resolve()

    sbatch_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --cpus-per-task={int(cpus)}",
        f"#SBATCH --mem={int(mem_gb)}G",
        f"#SBATCH --time={time}",
        f"#SBATCH --chdir={workdir}",
        f"#SBATCH --output={workdir}/{job_name}.%j.out",
        f"#SBATCH --error={workdir}/{job_name}.%j.err",
    ]
    if partition.strip():
        sbatch_lines.insert(2, f"#SBATCH --partition={partition.strip()}")

    # IMPORTANT: adjust if you need to load modules/conda here
    sbatch_lines += [
        "",
        "set -euo pipefail",
        "",
        "echo \"Running on $(hostname)\"",
        "echo \"Start: $(date)\"",
        " ".join(shlex.quote(c) for c in cmd_list),
        "echo \"End: $(date)\"",
    ]
    script_text = "\n".join(sbatch_lines)

    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".sbatch") as f:
        f.write(script_text)
        script_path = f.name

    st.code(script_text, language="bash")

    # submit
    proc = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
    if proc.returncode == 0:
        st.success(proc.stdout.strip())
        st.info(f"Sbatch script: {script_path}")
    else:
        st.error("sbatch failed")
        st.code(proc.stdout + "\n" + proc.stderr)

# ---------- Global options ----------
st.sidebar.header("Execution mode")
mode = st.sidebar.selectbox("How to run?", ["Local run (Popen)", "Slurm submit (sbatch)"])

st.sidebar.header("Common")
workdir = st.sidebar.text_input("Working directory (optional)", value="")
extra_global = st.sidebar.text_input("Extra args (appended verbatim)", value="")

# Slurm options
if mode == "Slurm submit (sbatch)":
    st.sidebar.header("Slurm resources")
    slurm_job = st.sidebar.text_input("Job name", value="snakehichip_ui")
    slurm_cpus = st.sidebar.number_input("CPUs", min_value=1, max_value=512, value=64)
    slurm_mem = st.sidebar.number_input("Memory (GB)", min_value=1, max_value=4096, value=128)
    slurm_time = st.sidebar.text_input("Time (HH:MM:SS)", value="24:00:00")
    slurm_partition = st.sidebar.text_input("Partition (optional)", value="")

def execute(cmd, job_name):
    if extra_global.strip():
        cmd += shlex.split(extra_global)

    if mode == "Local run (Popen)":
        run_and_stream_local(cmd, workdir=workdir or None)
    else:
        submit_slurm(
            cmd,
            job_name=job_name,
            cpus=slurm_cpus,
            mem_gb=slurm_mem,
            time=slurm_time,
            partition=slurm_partition,
            workdir=workdir or None,
        )

tabs = st.tabs(["Genomesetup", "Snakehichip", "Diffhichip", "TFscprinter", "TFtobias"])

# ---------- Tab 1: Genomesetup ----------
with tabs[0]:
    st.subheader("Genomesetup — build genome index/annotations")
    with st.form("genomesetup_form"):
        organism = st.text_input("--organism", value="hg38")
        fasta = st.text_input("--organism_fasta", value="/path/to/genomes/hg38/genome.fa")
        out_dir = st.text_input("-o", value="/path/to/output/genome")

        restriction_enzyme = st.text_input("--restriction_enzyme", value="mboi")
        bin_size = st.number_input("--bin_size", min_value=100, max_value=5000000, value=5000, step=1000)
        threads = st.number_input("--threads", min_value=1, max_value=512, value=150)

        submitted = st.form_submit_button("Start Genomesetup", type="primary")

    if submitted:
        cmd = [
            "Genomesetup",
            "--organism", organism,
            "--organism_fasta", fasta,
            "-o", out_dir,
            "--restriction_enzyme", restriction_enzyme,
            "--bin_size", str(int(bin_size)),
            "--threads", str(int(threads)),
        ]
        execute(cmd, job_name="genomesetup")

# ---------- Tab 2: Snakehichip ----------
with tabs[1]:
    st.subheader("Snakehichip — main pipeline")
    with st.form("snakehichip_form"):
        input_dir = st.text_input("-i", value="/path/to/input/HiCHIP_fastqs")
        genome_dir = st.text_input("--genome_dir", value="/path/to/output/genome")
        out_dir = st.text_input("-o", value="/path/to/output/HiCHIP_processed")

        organism = st.text_input("--organism", value="hg38")
        bin_size = st.number_input("--bin_size", min_value=100, max_value=5000000, value=5000, step=1000)
        threads = st.number_input("--threads", min_value=1, max_value=512, value=64)

        submitted = st.form_submit_button("Run Snakehichip", type="primary")

    if submitted:
        cmd = [
            "Snakehichip",
            "-i", input_dir,
            "--genome_dir", genome_dir,
            "-o", out_dir,
            "--threads", str(int(threads)),
            "--organism", organism,
            "--bin_size", str(int(bin_size)),
        ]
        execute(cmd, job_name="snakehichip")

# ---------- Tab 3: Diffhichip ----------
with tabs[2]:
    st.subheader("Diffhichip — differential analysis")
    with st.form("diffhichip_form"):
        sample_sheet = st.text_input("-s samplesheet", value="/path/to/samplesheet.tsv")
        out_dir = st.text_input("-o", value="/path/to/output/diffhichip")

        threads = st.number_input("--threads", min_value=1, max_value=512, value=2)
        bin_size = st.number_input("--bin_size", min_value=100, max_value=5000000, value=5000, step=1000)
        fitType = st.selectbox("--fitType", ["mean", "local", "parametric"], index=0)
        fdr = st.number_input("--fdr", min_value=0.0, max_value=1.0, value=0.05, step=0.01, format="%.3f")

        submitted = st.form_submit_button("Run Diffhichip", type="primary")

    if submitted:
        cmd = [
            "Diffhichip",
            "-s", sample_sheet,
            "-o", out_dir,
            "--threads", str(int(threads)),
            "--bin_size", str(int(bin_size)),
            "--fitType", fitType,
            "--fdr", str(fdr),
        ]
        execute(cmd, job_name="diffhichip")

# ---------- Tab 4: TFscprinter ----------
with tabs[3]:
    st.subheader("TFscprinter")
    with st.form("tfscprinter_form"):
        out_dir = st.text_input("-o", value="/path/to/output/tfscprinter")
        organism = st.text_input("--organism", value="hg38")
        atac_bam = st.text_input("--atacseq_bam", value="/path/to/input/atac.filtered.bam")
        enhancer_bedpe = st.text_input("--enhancer_bedpe", value="/path/to/input/significant_interactions.bedpe")
        samples = st.text_input("--samples", value="Sample1")
        threads = st.number_input("--threads", min_value=1, max_value=512, value=50)

        submitted = st.form_submit_button("Run TFscprinter", type="primary")

    if submitted:
        cmd = [
            "TFscprinter",
            "-o", out_dir,
            "--organism", organism,
            "--atacseq_bam", atac_bam,
            "--enhancer_bedpe", enhancer_bedpe,
            "--samples", samples,
            "--threads", str(int(threads)),
        ]
        execute(cmd, job_name="tfscprinter")

# ---------- Tab 5: TFtobias ----------
with tabs[4]:
    st.subheader("TFtobias")
    with st.form("tftobias_form"):
        out_dir = st.text_input("-o", value="/path/to/output/tftobias")
        genome_dir = st.text_input("--genome_dir", value="/path/to/output/genome")
        organism = st.text_input("--organism", value="hg38")
        atac_bam = st.text_input("--atacseq_bam", value="/path/to/input/atac.filtered.bam")
        enhancer_bed = st.text_input("--enhancer_bed", value="/path/to/input/significant_interactions.bedpe")
        samples = st.text_input("--samples", value="Sample1")
        motif = st.text_input("--motif", value="/path/to/motifs/JASPAR2024_CORE_non-redundant_pfms.meme")
        threads = st.number_input("--threads", min_value=1, max_value=512, value=50)

        submitted = st.form_submit_button("Run TFtobias", type="primary")

    if submitted:
        cmd = [
            "TFtobias",
            "-o", out_dir,
            "--genome_dir", genome_dir,
            "--organism", organism,
            "--atacseq_bam", atac_bam,
            "--enhancer_bed", enhancer_bed,
            "--samples", samples,
            "--threads", str(int(threads)),
            "--motif", motif,
        ]
        execute(cmd, job_name="tftobias")