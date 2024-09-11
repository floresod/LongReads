import glob
import os 

SAMPLE,=glob_wildcards("../resources/Data/Fastq/{sample}.fastq")

rule all:
    input:
        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE), 
        "../results/multiqc_rawreads/multiqc_report.html"

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

