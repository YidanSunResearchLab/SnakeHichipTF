import os

#############################################
# Configuration
#############################################
# Configuration
script_dir = os.path.dirname(workflow.snakefile)
output_dir = config["output_dir"]
organism = config["organism"]
organism_fasta = config["organism_fasta"]
restriction_enzyme = config["restriction_enzyme"]
bin_size = config["bin_size"]
threads = int(config["threads"] / 2)
read_length = 50  # For mappability
frag_size = 500


# Embedded config for restriction enzymes
config["re"] = {
    "hindiii": {"site": "AAGCTT", "pos": "1", "hicpro": "A^AGCTT", "ligation": "AAGCTTAAGCTT"},
    "mboi": {"site": "GATC", "pos": "0", "hicpro": "^GATC", "ligation": "GATCGATC"},
    "ncoi": {"site": "CCATGG", "pos": "1", "hicpro": "C^CATGG", "ligation": "CCATGGCCATGG"},
    "arima": {"site": "GATC,GANTC", "pos": "0,1", "hicpro": "^GATC,G^ANTC", "ligation": "GATCGATC,GATCGANT,GANTGATC,GANTGANT"}, #GATCGATC,GATCGANT,GANTGATC,GANTGANT 
    "mnase": {"site": "mnase", "pos": "none", "hicpro": "mnase", "ligation": None},
    "dpnii": {"site": "GATC", "pos": "0", "hicpro": "^GATC", "ligation": "GATCGATC"},
    "bglii": {"site": "AGATCT", "pos": "1", "hicpro": "A^GATCT", "ligation": "AGATCTAGATCT"},
    "msei": {"site": "TTAA", "pos": "1", "hicpro": "T^TAA", "ligation": "TTAATTAA"},
    "hinfI": {"site": "GANTC", "pos": "1", "hicpro": "G^ANTC", "ligation": "GANTCGANTC"},
    "nlaIII": {"site": "CATG", "pos": "4", "hicpro": "CATG^", "ligation": "CATGCATG"}
}

# Check if genome or restriction enzyme setup is needed
GENOME_SETUP_NEEDED = not os.path.exists(output_dir + "/organisms/" + organism + "/" + organism + ".txt") and organism_fasta != "NA"
MAPS_SETUP_NEEDED = not os.path.exists(output_dir + "/organisms/" + organism + "/" +restriction_enzyme + "/" + str(bin_size) + "/restriction_enzyme.organism.txt")
HICDC_SETUP_NEEDED = not os.path.exists(output_dir + "/organisms/" + organism + "/" +restriction_enzyme + "/" + str(bin_size) + "/restriction_enzyme.hicdcplus.bintolen.txt.gz")
Fithichip_SETUP_NEEDED = not os.path.exists(output_dir + "/organisms/" + organism + "/" +restriction_enzyme + "/" + str(bin_size) + "/restriction_enzyme.fithichip.bed")



#############################################
# Final target
#############################################
rule all:
    input:
        expand(output_dir + "/organisms/{organism}/config_used.yaml", organism=organism),
        expand(output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt", organism=organism) if GENOME_SETUP_NEEDED else [], #Genome Index Generation
        expand(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.organism.txt", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size) if MAPS_SETUP_NEEDED else [], #Genomic Features Generation
        expand(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.hicdcplus.bintolen.txt.gz", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size) if HICDC_SETUP_NEEDED else [],
        expand(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.fithichip.bed", organism=organism, restriction_enzyme=restriction_enzyme, bin_size=bin_size) if Fithichip_SETUP_NEEDED else [],


#############################################
# Config file output
#############################################
rule show_config:
    output: output_dir + "/organisms/{organism}/config_used.yaml"
    run:
        import yaml
        with open(output[0], "w") as f:
            yaml.dump(config, f, sort_keys=True)

#############################################
# --- Genome Index Generation Rules ---
#############################################
#if not os.path.exists(output_dir + "/organisms/" + organism + "/" + organism + ".txt") and organism_fasta != "NA":
def get_main_chromosomes(organism, fai_file):
    """Infer main chromosomes based on organism and FASTA index."""
    # Predefined patterns for common organisms (adjust as needed)
    chrom_patterns = {
        "hg38": r"^chr([1-9][0-9]*|[XY])$",  # chr1-22, chrX, chrY
        "mm10": r"^chr([1-9][0-9]*|[XY])$",  # chr1-19, chrX, chrY
        "hg19": r"^chr([1-9][0-9]*|[XY])$",  # chr1-22, chrX, chrY
        "mm9": r"^chr([1-9][0-9]*|[XY])$",  # chr1-19, chrX, chrY
        "dm6": r"^chr(2L|2R|3L|3R|4|X|Y)$",  # Drosophila main arms
        "dm3": r"^chr(2L|2R|3L|3R|4|X|Y)$"  # Drosophila main arms
    }
    pattern = chrom_patterns.get(organism, r"^chr([1-9][0-9]*|[XY])$")  # Default to mammal pattern
    import re
    main_chroms = []
    if os.path.exists(fai_file):
        with open(fai_file, "r") as f:
            for line in f:
                chrom = line.split("\t")[0]
                if re.match(pattern, chrom):
                    main_chroms.append(chrom)
    return main_chroms if main_chroms else ["chr1"]  # Fallback to avoid empty list

# Rule to copy FASTA file and create organism directory
rule copy_fasta:
    output:
        fasta=output_dir + "/organisms/{organism}/{organism}.fa",
        fai=output_dir + "/organisms/{organism}/{organism}.fa.fai"  # samtools index 
    params:
        organism_dir=output_dir + "/organisms/{organism}",
        fasta=organism_fasta if organism_fasta!="NA" else []
    threads: threads
    conda: "envs/samtools.yaml"
    log:
        "logs/copy_fasta_{organism}.log"
    benchmark:
        "benchmarks/copy_fasta_{organism}.bmk.txt"
    shell:
        """
        if [ "{params.fasta}" == "" ]; then
            echo "Error: organism_fasta is 'NA' and no FASTA file provided." >&2
            exit 1
        fi
        mkdir -p {params.organism_dir}
        cp {params.fasta} {output.fasta}
        samtools faidx {output.fasta}
        python -c "with open('{output.fai}') as f: \
                main_chroms = [line.split()[0] for line in f if line.split()[0].startswith('chr')]; \
                print(' '.join(main_chroms))" > {params.organism_dir}/main_chroms.txt 2>> {log}
        samtools faidx {output.fasta} $(cat {params.organism_dir}/main_chroms.txt) > {params.organism_dir}/{wildcards.organism}_filtered.fa
        mv {params.organism_dir}/{wildcards.organism}_filtered.fa {output.fasta}
        samtools faidx {output.fasta}  # Reindex after filtering
        """

# Rule to get chromosome size
rule generate_chrom_sizes:
    input:
        fai=output_dir + "/organisms/{organism}/{organism}.fa.fai"  # Temp index for filtering
    output:
        output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    conda:
        "envs/basic.yaml"
    log:
        "logs/generate_chrom_sizes/{organism}.log"
    benchmark:
        "benchmarks/generate_chrom_sizes/{organism}.txt"
    shell:
        """
        cut -f 1,2 {input.fai} > {output} 2> {log}
        """

# Rule to calculate autosome count
rule calculate_autosome_count:
    input:
        chrom_sizes = output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    output:
        output_dir+"/organisms/{organism}/{organism}.autosome.count.txt"
    run:
        with open(input.chrom_sizes, "r") as f:
            autosomes = sum(1 for line in f if line.split()[0].startswith("chr") and 
                            line.split()[0] not in ["chrX", "chrY", "chrM"])
        with open(output[0], "w") as out:
            out.write(str(autosomes))

# Rule to build Bowtie2 index
rule bowtie2Index:
    input:
        output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    output:
        output_dir+"/organisms/{organism}/BowtieIndex/{organism}.rev.2.bt2"
    params:
        basedir=output_dir+"/organisms/{organism}/BowtieIndex",
        genome = output_dir+"/organisms/{organism}/{organism}.fa"
    threads: threads
    conda: "envs/basic.yaml"
    log:
        "logs/bowtie2Index_{organism}.log"
    benchmark:
        "benchmarks/bowtie2Index_{organism}.bmk.txt"
    shell:
        """
        mkdir -p {params.basedir}
        bowtie2-build --threads {threads} {params.genome} {params.basedir}/{organism} > {log} 2>&1
        """

# Rule to build Chromap index
rule chromapIndex:
    input:
        output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    output:
        output_dir + "/organisms/{organism}/ChromapIndex/{organism}.index"
    params:
        genome = output_dir + "/organisms/{organism}/{organism}.fa",  # Input FASTA
        index_dir = output_dir + "/organisms/{organism}/ChromapIndex",
        index_file = output_dir + "/organisms/{organism}/ChromapIndex/{organism}.index"
    threads: threads
    conda: "envs/maps.yaml"  # Assuming chromap is available here
    log:
        "logs/chromapIndex_{organism}.log"
    benchmark:
        "benchmarks/chromapIndex_{organism}.bmk.txt"
    shell:
        """
        mkdir -p {params.index_dir}
        chromap -i -r {params.genome} -o {params.index_file} -t {threads} > {log} 2>&1
        """
        
# Calculate mappability with genmap
# Rule to index the genome with genmap index
rule genmap_index:
    input:
        output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    output:
        index=output_dir+"/organisms/{organism}/GenmapIndex.done",
    params:
        index_dir=output_dir+"/organisms/{organism}/GenmapIndex/",  # Directory for index files
        genome=output_dir + "/organisms/{organism}/{organism}.fa"
    conda:
        "envs/basic.yaml"
    log:
        "logs/GenmapIndex.{organism}.log"
    benchmark:
        "benchmarks/GenmapIndex.{organism}.txt"
    threads: threads
    shell:
        """
        if [ -d "{params.index_dir}" ]; then
            rm -rf {params.index_dir} 
        fi
        genmap index -F {params.genome} -I {params.index_dir} > {log} 2>&1 && 
        touch {output}
        """

# Rule to compute mappability with genmap map
rule genmap_mappability:
    input:
        index=output_dir+"/organisms/{organism}/GenmapIndex.done",
    output:
        wig=output_dir+"/organisms/{organism}/{organism}_mappability_" + str(read_length) + ".wig"
    params:
        index_dir=output_dir+"/organisms/{organism}/GenmapIndex/",  # Directory for index files
        output_prefix=output_dir+"/organisms/{organism}/{organism}_mappability_" + str(read_length) + "",
    conda:
        "envs/basic.yaml"
    log:
        "logs/genmap_mappability.{organism}_" + str(read_length) + ".log"
    benchmark:
        "benchmarks/genmap_mappability.{organism}_" + str(read_length) + ".txt"
    threads: threads
    shell:
        """
        genmap map -I {params.index_dir} -K {read_length} -E 2 -w -O {params.output_prefix} -T {threads} > {log} 2>&1
        """

# Rule to convert WIG to BigWig
rule wig_to_bigwig:
    input:
        wig=output_dir+"/organisms/{organism}/{organism}_mappability_" + str(read_length) + ".wig",
        sizes=output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    output:
        bigwig=output_dir+"/organisms/{organism}/{organism}_mappability_" + str(read_length) + ".bw"
    params:
        wig=output_dir+"/organisms/{organism}/{organism}_mappability_" + str(read_length) + ".wig",
        sizes=output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    conda:
        "envs/genomefeature.yaml"
    log:
        "logs/wig_to_bigwig.{organism}_" + str(read_length) + ".log"
    benchmark:
        "benchmarks/wig_to_bigwig.{organism}_" + str(read_length) + ".txt"
    threads: threads
    shell:
        """
        wigToBigWig {params.wig} {params.sizes} {output.bigwig} > {log} 2>&1
        """

# Rule to modify basic.txt
rule modify_basic_txt:
    input:
        genome=output_dir+"/organisms/{organism}/{organism}.fa",
        fa_index=output_dir+"/organisms/{organism}/{organism}.fa.fai",
        chrom_size=output_dir+"/organisms/{organism}/{organism}.chrom.sizes",
        autosome=output_dir+"/organisms/{organism}/{organism}.autosome.count.txt",
        bowtie2_index=output_dir+"/organisms/{organism}/BowtieIndex/{organism}.rev.2.bt2",
        chromap_index=output_dir + "/organisms/{organism}/ChromapIndex/{organism}.index",
        bigwig=output_dir+"/organisms/{organism}/{organism}_mappability_" + str(read_length) + ".bw"
    output:
        output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt"
    params:
        basic_txt = script_dir+"/../organisms/basic.txt",
        ligation_site = config["re"][restriction_enzyme]["ligation"],
    conda: "envs/basic.yaml"
    log:
        "logs/modify_basic_txt_{organism}.log"
    benchmark:
        "benchmarks/modify_basic_txt_{organism}.bmk.txt"
    threads: min(threads, 64)
    shell:
        """
        sed 's|GENOME_FRAGMENT = organism_dir|GENOME_FRAGMENT = {output_dir}/organisms/{wildcards.organism}/{restriction_enzyme}/{bin_size}|g' {params.basic_txt} \
        | sed 's|organism_dir|{output_dir}/organisms/{wildcards.organism}|g' \
        | sed 's|threads_input|{threads}|g' \
        | sed 's|genomename|{wildcards.organism}|g' \
        | sed 's|BIN_SIZE =|BIN_SIZE = {bin_size}|g' \
        | sed 's|LIGATION_SITE =|LIGATION_SITE = {params.ligation_site}|g' \
        > {output}.tmp 2> {log}

        # Add --trim5 5 for Arima only
        if [ "{restriction_enzyme}" = "arima" ]; then
            sed 's|--very-sensitive|--trim5 5 --very-sensitive|g' {output}.tmp > {output}
        else
            mv {output}.tmp {output}
        fi
        """


#############################################
# --- Restriction Enzyme Fragments Generation Rules ---
#############################################
rule find_resfrags:
    input:
        annotation=output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt"
    output:
        output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_resolution_fragments.bed"
    params:
        re_site=lambda wildcards: config["re"][wildcards.restriction_enzyme]["hicpro"],
        genome=output_dir+"/organisms/{organism}/{organism}.fa",
    log:
        "logs/find_resfrags.{organism}.{restriction_enzyme}.{bin_size}.restriction_enzyme_resolution_fragments.log"
    benchmark:
        "benchmarks/find_resfrags.{organism}.{restriction_enzyme}.{bin_size}.restriction_enzyme_resolution_fragments.txt"
    conda:
        "envs/hicpro.yaml"  # Conda environment with HiC-Pro installed
    shell:
        """
        {script_dir}/../scripts/hicpro-3.1.0/bin/utils/digest_genome.py -r {params.re_site} -o {output} {params.genome} > {log} 2>&1
        """

# Additional rule for RE_cut_multiple (from MAPS), this step identify restriction enzyme (RE) cut sites or, alternatively, to segment the genome into fixed-size bins (for MNase-based analysis)
#RE Mode:
#Identifies all occurrences of the RE sequence(s) in the genome.
#Cuts at specified positions, creating fragments between consecutive sites.
#For each cut site, generates two fragments (one on the negative strand, one on the positive strand) with their end/start positions, lengths, and GC content.
#MNase Mode:
#Divides the genome into fixed-size bins (e.g., 5000 bp).
#Assigns fragment-like entries for each bin, mimicking RE fragments.
rule RE_cut_multiple:
    input:
        annotation=output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt"
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.txt"
    params:
        genome=output_dir+"/organisms/{organism}/{organism}.fa",
        cut=script_dir+"/../scripts/restriction_cut_multipleenzyme.py",
        site=lambda wildcards: config["re"][restriction_enzyme]["site"],
        pos=lambda wildcards: config["re"][restriction_enzyme]["pos"],
        bin_size=bin_size
    conda:
        "envs/genomefeature.yaml"
    log:
        "logs/re_cut_multiple.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/re_cut_multiple.{organism}.{restriction_enzyme}.{bin_size}.txt"
    shell:
        """
        python {params.cut} -f {params.genome} -s {params.site} -p {params.pos} -o {output} -b {params.bin_size} > {log} 2>&1
        """

# Get fragment ends (fends)
rule get_fend:
    input:
        cut=output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.txt",
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.bed"
    params:
        frag_size=frag_size,
        genome_size=output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    conda:
        "envs/basic.yaml"
    log:
        "logs/get_fend.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/get_fend.{organism}.{restriction_enzyme}.{bin_size}.txt"
    shell:
        """
        awk -vOFS="\\t" '{{print $3,$4,$4,$3"_"$1,"0",$2}}' {input.cut} | \
        bedtools slop -s -l 0 -r {params.frag_size} -g {params.genome_size} > {output} 2> {log}
        """

rule cal_mappability:
    input:
        annotation=output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt",
        fend=output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.bed"
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.map.tab"
    params:
        mapability=output_dir+"/organisms/{organism}/"+organism+"_mappability_" + str(read_length) + ".bw",
    conda:
        "envs/genomefeature.yaml"
    log:
        "logs/cal_mappability.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/cal_mappability.{organism}.{restriction_enzyme}.{bin_size}.txt"
    shell:
        """
        bigWigAverageOverBed {params.mapability} {input.fend} {output} 2> {log}
        """

# Merge mappability
rule merge_map:
    input:
        map=output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.map.tab",
        cut=output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.txt"
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.map"
    params:
        merge=script_dir+"/../scripts/merge_map.py"
    conda:
        "envs/genomefeature.yaml"
    log:
        "logs/merge_map.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/merge_map.{organism}.{restriction_enzyme}.{bin_size}.txt"
    shell:
        """
        python {params.merge} -c {input.cut} -m {input.map} -o {output} > {log} 2>&1
        """

# Filter mappability
rule filter_map:
    input:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.map"
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.map.filter"
    conda:
        "envs/genomefeature.yaml"
    log:
        "logs/filter_map.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/filter_map.{organism}.{restriction_enzyme}.{bin_size}.txt"
    shell:
        """
        awk '$7>0.5 && $5 <= 20000' {input} > {output} 2> {log}
        """

rule bin_feature:
    input:
        feat=output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_position_all.map.filter",
        hicpro=output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_resolution_fragments.bed"
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.organism.txt"
    params:
        genome_size=expand(output_dir+"/organisms/{organism}/{organism}.chrom.sizes", organism=organism),
        feat2bin=script_dir+"/../scripts/feature_frag2bin.py",
        bin_size=bin_size
    conda:
        "envs/genomefeature.yaml"
    log:
        "logs/bin_feature.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/bin_feature.{organism}.{restriction_enzyme}.{bin_size}.txt"
    shell:
        """
        python {params.feat2bin} -i {input.feat} -o {output} -b {params.bin_size} -g {params.genome_size} > {log} 2>&1
        """


# Rule to generate genomic features for hicdcplus
rule hicdc_genomefeature:
    input:
        annotation=output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt"
    output:
        output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.hicdcplus.bintolen.txt.gz"
    params:
        outdir=output_dir+"/organisms/{organism}/{restriction_enzyme}/{bin_size}",
        enzyme=lambda wildcards: config["re"][restriction_enzyme]["site"],
        genome=output_dir + "/organisms/" + organism + "/" + organism + ".fa",
        binsize=bin_size,
        chrom_sizes=output_dir+"/organisms/{organism}/{organism}.chrom.sizes"
    log:
        "logs/hicdc_genomefeature.{organism}.{restriction_enzyme}.{bin_size}.log"
    benchmark:
        "benchmarks/hicdc_genomefeature.{organism}.{restriction_enzyme}.{bin_size}.txt"
    threads: threads
    conda:
        "envs/hicdc.yaml"  
    shell:
        """
        Rscript {script_dir}/../scripts/hicdc_genomefeature.R \
          -o {params.outdir} \
          -g {params.genome} \
          -r {params.enzyme} \
          -b {params.binsize} \
          -c {params.chrom_sizes} > {log} 2>&1
        """


# Rule to generate genomic features for fithichip
rule fithichip_mappability:
    input:
        annotation=output_dir + "/organisms/{organism}/"+restriction_enzyme+"/"+str(bin_size)+"/{organism}.txt",
    output:
        mappability = output_dir+"/organisms/{organism}/"+organism+"_mappability_" + str(read_length) + ".bedgraph",
    params:
        mappability = output_dir+"/organisms/{organism}/"+organism+"_mappability_" + str(read_length) + ".bw",
    log:
        "logs/fithichip_mappability.{organism}.log"
    conda:
        "envs/genomefeature.yaml"  
    shell:
        """
        bigWigToBedGraph {params.mappability} {output.mappability}
        """

rule fithichip_features:
    input:
        re_frags = output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme_resolution_fragments.bed",
        mappability = output_dir+"/organisms/{organism}/"+organism+"_mappability_" + str(read_length) + ".bedgraph",
    output:
        features = output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/restriction_enzyme.fithichip.bed",
        temp_mapp_offset = temp(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/Temp_Fragment_Mapp_200bp.bed"),
        temp_mapp = temp(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/Mappability_RE_Fragments.bed"),
        temp_gc_offset = temp(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/Temp_Fragment_GC_200bp.bed"),
        temp_gc = temp(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/GC_Content_RE_Fragments.bed"),
        temp_mapp_dump = temp(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/Mappability_Dump.bed"),
        temp_gc_dump = temp(output_dir + "/organisms/{organism}/{restriction_enzyme}/{bin_size}/GC_Dump.bed")
    params:
        ref_fasta = output_dir + "/organisms/{organism}/{organism}.fa", 
        mapp_window = 200, 
        gc_window = 200,    
    log:
        "logs/fithichip_features.{organism}.{restriction_enzyme}.{bin_size}.log"
    conda:
        "envs/fithichip.yaml"
    shell:
        """
        echo "======== Creating the fragment end (w.r.t window size) file -- to compute the mappability information!!" > {log}
        awk -v s={params.mapp_window} -F'\\t' 'function max(x,y){{return x>y?x:y}}; function min(x,y){{return x<y?x:y}}; {{printf "%s\\t%d\\t%d\\n%s\\t%d\\t%d\\n", $1, $2, min($2+s,$3), $1, max($3-s,$2), $3}}' {input.re_frags} | sort -k1,1 -k2,2n > {output.temp_mapp_offset} 2>> {log}

        echo "======== Creating the mappability file !!" >> {log}
        bedtools map -a {output.temp_mapp_offset} -b {input.mappability} -c 4 -o mean | awk -F'\\t' '{{if ($4=="." || $4=="NA" || $4=="NaN") {{$4=0}}; print $0}}' > {output.temp_mapp} 2>> {log}

        echo "======== Creating the fragment end (w.r.t window size) file -- to compute the GC content information!!" >> {log}
        awk -v s={params.gc_window} -F'\\t' 'function max(x,y){{return x>y?x:y}}; function min(x,y){{return x<y?x:y}}; {{printf "%s\\t%d\\t%d\\n%s\\t%d\\t%d\\n", $1, $2, min($2+s,$3), $1, max($3-s,$2), $3}}' {input.re_frags} | sort -k1,1 -k2,2n > {output.temp_gc_offset} 2>> {log}

        echo "======== Creating the GC content file !!" >> {log}
        bedtools nuc -fi {params.ref_fasta} -bed {output.temp_gc_offset} | awk -F'\\t' '{{if ($4=="." || $4=="NA" || $4=="NaN") {{$4=0}}; print $0}}' > {output.temp_gc} 2>> {log}

        awk -F'\\t' '{{{{sum+=$4}}}} NR%2==0 {{{{print sum/2; sum=0}}}}' {output.temp_mapp} > {output.temp_mapp_dump} 2>> {log}
        awk -F'\\t' '{{if (NR>1) {{{{sum+=$5}}}}}}; {{{{if (NR%2!=0 && NR>1) {{{{print sum/2; sum=0}}}}}}}}' {output.temp_gc} > {output.temp_gc_dump} 2>> {log}

        Rscript {script_dir}/../scripts/fithichip/CombineREFragMappGC.r {input.re_frags} {output.temp_mapp_dump} {output.temp_gc_dump} {output.features} >> {log} 2>&1
        """

