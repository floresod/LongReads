import os
import glob 

SAMPLE,=glob_wildcards("../resources/Data/Fastq/{sample}.fastq")

rule all: 
    input:
        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE),
        "../results/multiqc_rawreads/multiqc_report.html"
#        expand("Kraken_report/{sample}.txt", sample=SAMPLE), 
#        expand("Bracken_report/{sample}.txt", sample=SAMPLE),
#        expand("Contigs/flye/{sample}/assembly.fasta", sample=SAMPLE), 
#        expand("Kraken_contigs/{sample}.txt", sample=SAMPLE), 
#        expand("Kraken_contigs_out/{sample}.txt", sample=SAMPLE)

rule fastqc_rawreads: 
    input: 
        rawread="../resources/Data/Fastq/{sample}.fastq"
    output: 
        html="../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html"
    conda:
        "../envs/fastqc_env.yaml"
    log:
        "../resources/Logs/fastqc_rawreads/{sample}.log"
    threads:
        4
    params:
        path="../resources/Outputs/fastqc_rawreads/",
        memory="4000"
    shell:
        """
        fastqc  {input.rawread} \
                --threads {threads} \
                --memory {params.memory} \
                -o {params.path}  > {log} 2>&1
        """ 

rule multiqc_rawreads:
    input: 
        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE)
    output: 
       "results/multiqc_rawreads/multiqc_report.html", 
        directory("../results/multiqc_data")
       log: 
        "../resources/Logs/multiqc_rawreads/multiqc_rawreads.log"
#    conda:
#       "../envs/fastqc_env.yaml"
#   shell: 
#        """
#        mkdir -p {params.path}
#
#        multiqc {input} -o {params.path} \
#        --filename {output.report} > {log} 2>&1
#        """
    wrapper:
        "v4.3.0/bio/multiqc"

rule kraken2_run: 
    input: 
        rawread="../resources/Data/Fastq/{sample}.fastq"
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
                {input.rawread} \
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

        

