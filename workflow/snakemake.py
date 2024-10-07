import glob
import os 

SAMPLE,=glob_wildcards("../resources/Data/Fastq/{sample}.fastq")

rule all:
    input:
        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE), 
        expand("../results/kraken2_rr/{sample}.txt", sample=SAMPLE),
        expand("../results/bracken/{sample}.txt", sample=SAMPLE),
        expand("../resources/Outputs/flye/{sample}/assembly.fasta", sample=SAMPLE),
        expand("../resources/Outputs/medaka/{sample}/consensus.fasta", sample=SAMPLE), 
        expand("../results/FinalContigs/{sample}.fasta", sample = SAMPLE),
        expand("../results/checkm2/{sample}.tsv", sample=SAMPLE)

rule fastqc_rawreads:
    input:
        "../resources/Data/Fastq/{sample}.fastq"
    output:
        html="../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html",
        zip="../resources/Outputs/fastqc_rawreads/{sample}_fastqc.zip"
    params:
        exra="--quiet"
    log:
        "../resources/Logs/fastqc_rawreads/{sample}.log"
    threads:
        4
    resources:
        mem_mb=4000
    wrapper:
        "v4.6.0/bio/fastqc"

rule kraken2_rr:
    input:
        rawread="../resources/Data/Fastq/{sample}.fastq"
    output:
        report="../results/kraken2_rr/{sample}.txt", 
        output="../resources/Outputs/kraken2_rr/{sample}.txt"
    params:
        database="../../../Databases/k2_standard_08gb_20240605", 
        threads="6"
    log:
        "../resources/Logs/kraken2_rawreads/{sample}.log"
    conda: 
        "../envs/kraken2_env.yaml"
    shell:
        """
        kraken2 --db {params.database} \
                {input.rawread} \
                --report {output.report} \
                --output {output.output} \
                --threads {threads} \
                --confidence 0.005 \
                --use-names > {log} 2>&1
        """

rule bracken_rr:
    input: 
        report="../results/kraken2_rr/{sample}.txt"
    output:
        report="../results/bracken/{sample}.txt"
    params:
        database="../../../Databases/k2_standard_08gb_20240605", 
        length=100
    log: 
        "../resources/Logs/bracken/{sample}.log"
    conda:
        "../envs/bracken_env.yaml"
    shell:
        """
        bracken -d {params.database} \
                -i {input.report} \
                -o {output.report} \
                -r {params.length} > {log} 2>&1
        """

rule assembly_flye:
    input: 
        reads = "../resources/Data/Fastq/{sample}.fastq"
    output:
        contigs = "../resources/Outputs/flye/{sample}/assembly.fasta" 
    params:
        outdir = "../resources/Outputs/flye/{sample}",
        genome_size = "2g"
    threads:
        6
    conda:
        "../envs/flye_env.yaml"
    log: 
        "../resources/Logs/flye/{sample}.log"
    shell: 
        """
        mkdir -p {params.outdir}
        flye    --nano-raw {input.reads} \
                --out-dir {params.outdir} \
                --genome-size {params.genome_size}\
                --threads {threads}\
                --meta  > {log} 2>&1
        """

rule polish_medaka:
    input:
        assembly="../resources/Outputs/flye/{sample}/assembly.fasta",
        basecall="../resources/Data/Fastq/{sample}.fastq"
    output:
        "../resources/Outputs/medaka/{sample}/consensus.fasta"
    params:
        outdir="../resources/Outputs/medaka/{sample}"
    threads: 
        6
    log:
        "../resources/Logs/medaka/{sample}.log"
    conda:
        "../envs/medaka_env.yaml"
    shell:
        """
        medaka_consensus -i {input.basecall} \
                         -d {input.assembly} \
                         -o {params.outdir} \
                         -t {threads} > {log} 2>&1
        """
    
rule gathering_contigs:
    input:
        "../resources/Outputs/medaka/{sample}/consensus.fasta"
    output:
        "../results/FinalContigs/{sample}.fasta"
    log:
        "../resources/Logs/move_contigs/{sample}.log"
    shell: 
        """
        cp {input} {output}
        """
rule checkm2: 
    input: 
       contigs="../results/FinalContigs/{sample}.fasta"
    output:
        dir_report=directory("../resources/Outputs/checkm2/{sample}"),
        report="../results/checkm2/{sample}.tsv"
    log:
        "../resources/Logs/checkm2/{sample}.log"
    resources:
        threads=4
    conda:
        "../envs/checkm2_env.yaml"
    shell:
        """
        rm -rf {output.dir_report}

        checkm2 predict --input {input.contigs} \
                        --output-directory {output.dir_report} \
                        --force \
                        --threads {resources.threads} > {log} 2>&1
        
        mkdir -p "../results/checkm2"
        cp ../resources/Outputs/checkm2/{wildcards.sample}/quality_report.tsv \
            {output.report}
        """

        

