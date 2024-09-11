import glob
import os 

SAMPLE,=glob_wildcards("../resources/Data/Fastq/{sample}.fastq")

rule all:
    input:
        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE), 
#        "../results/multiqc_rawreads/multiqc_report.html"
        expand("../results/kraken2_rr/{sample}.txt", sample=SAMPLE),
        expand("../results/bracken/{sample}.txt", sample=SAMPLE)

rule fastqc_rawreads:
    input:
        rawread="../resources/Data/Fastq/{sample}.fastq"
    output:
        html="../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html"
    params:
        path="../resources/Outputs/fastqc_rawreads/", 
        memory=4000
    log:
        "../resources/Logs/fastqc_rawreads/{sample}.log"
    threads:
        4
    conda:
        "../envs/fastqc_env.yaml"
    shell:
        """
        fastqc  {input.rawread} \
                --threads {threads} \
                --memory {params.memory} \
                -o {params.path} > {log} 2>&1
        """

#rule multiqc_rawreads:
#    input: 
#        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE)
#    output:
#        directory("../results/multiqc_rawreads/"),
#        report="../results/multiqc_rawreads/multiqc_report.html"
#    params:
#        path="../results/multiqc_rawreads/"
#    log:    
#        "../resources/Logs/multiqc_rawreads/multiqc_rawreads.log"
#    conda:
#        "../envs/fastqc_env.yaml"
#    shell:
#        """
#        mkdir -p ../results/multiqc_rawreads
#        multiqc {input} -o {params.path} --filename {output.report}  > {log} 2>&1
#        """

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







