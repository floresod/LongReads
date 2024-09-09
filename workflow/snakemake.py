import os
import glob 

SAMPLE,=glob_wildcards("../resources/Fastq/{sample}.fastq.gz")

rule all: 
    input:
        "../results/multiqc_report_RawReads.html",
        expand()
#        expand("Kraken_report/{sample}.txt", sample=SAMPLE), 
#        expand("Bracken_report/{sample}.txt", sample=SAMPLE),
#        expand("Contigs/flye/{sample}/assembly.fasta", sample=SAMPLE), 
#        expand("Kraken_contigs/{sample}.txt", sample=SAMPLE), 
#        expand("Kraken_contigs_out/{sample}.txt", sample=SAMPLE)

rule fastqc_rawreads: 
    input: 
        rawread="../resources/Data/Fastq/{sample}.fastq.gz"
    output: 
        zip="../resources/Outputs/fastqc_rawreads/{sample}.zip",
        html="../resources/Outputs/fastqc_rawreads/{sample}.html"
    conda:
        "../envs/fastqc_env.yaml"
    log: 
        "../resources/Logs/fastqc_rawreads/{sample}.log"
    threads:
        4
    params:
        path="../resources/Outputs/fastqc_rawreads/"
    shell:
        """
        mkdir -p {params.path}
        mkdir -p ../resources/Logs/fastqc_rawreads

        fastqc  {input.rawread} \
                --threads {threads} \
                -o {params.path}  > {log} 2>&1
        """ 

rule multiqc_rawreads:
    input: 
        expand("../resources/Outputs/fastqc_rawreads/{sample}.html", sample=SAMPLE)

    output: 
        "../results/multiqc_report_RawReads.html" 
    log: 
        "../resources/Logs/multiqc_rawreads/multiqc_rawreads.log"
    conda:
        "../envs/fastqc_env.yaml"
    shell:
        """
        mkdir -p ../resources/Logs/multiqc_rawreads

        multiqc {input} -o {output} > {log} 2>&1
        """

rule kraken2_run: 
    input: 
        rawread="../resources/Data/Fastq/{sample}.fastq.gz"
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        threads=5
    output:
        output="../resources/Outputs/Kraken_output/{sample}.txt",
        report="../resources/Outputs/Kraken_report/{sample}.txt" 
    conda:
        "../envs/kraken2_env.yaml"
    log:
        "../resources/Logs/kraken2_{sample}.log"
    shell:
        """
        kraken2 --db {params.database} \
                --gzip-compressed {input.rawread} \
                --output {output.output} \
                --report {output.report} \
                --threads {params.threads} \
                --confidence 0.005 \
                --use-names   > {log} 2>&1 
        """

rule bracken_run:
    input:
        report=rules.kraken2_run.output.report
    output:
        report="../resources/Output/Bracken_report/{sample}.txt"
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        length=100
    conda:
        "../envs/bracken_env.yaml"
    log:
        "../resources/Logs/bracken_{sample}.log"

    shell:
        """
        bracken -d {params.database}\
                -i {input.report}\
                -o {output.report}\
                -r {params.length}  > {log} 2>&1
        """

rule assembly_flye:
    input: 
        reads = "../resources/Data/Fastq/{sample}.fastq.gz"
    output:
        contigs = "../resources/Outputs/flye/{sample}/assembly.fasta" 
    params:
        outdir = "../resources/Outputs/flye/{sample}",
        genome_size = "2g"
    threads:
        5
    conda:
        "../envs/flye_env.yaml"
    log: 
        "../resources/Logs/flye_{sample}.log"
    shell: 
        """
        mkdir -p {params.outdir}
        flye    --nano-raw {input.reads} \
                --out-dir {params.outdir} \
                --genome-size {params.genome_size}\
                --threads {threads}\
                --meta  > {log} 2>&1
        """

rule Kraken2_contigs: 
    input: 
        contigs = rules.assembly_flye.output
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        threads=5
    output:
        report="../results/Kraken_contigs/{sample}.txt",
        output="../results/Kraken_contigs_out/{sample}.txt"
    conda:
        "../envs/kraken2_env.yaml"
    log:
        "../resources/Logs/kraken2_Contigs_{sample}.log"
    shell:
        """
        kraken2 --db {params.database} {input.contigs} \
                --report {output.report} \
                --output {output.output} \
                --threads {params.threads} \
                --use-mpa-style \
                --use-names   > {log} 2>&1
        """

        

