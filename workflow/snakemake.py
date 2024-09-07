import os
import glob 

SAMPLE,=glob_wildcards("RawReads/{sample}.fastq.gz")

rule all: 
    input: 
        expand("RawQC/{sample}_fastqc.{extension}", sample=SAMPLE, extension=["html", "zip"]), 
#        expand("FilteredReads/{sample}.fastq.gz", sample=SAMPLE),
#        expand("FilteredQC/{sample}_fastqc.{extension}", sample=SAMPLE, extension=["html","zip"]),
        expand("Kraken_report/{sample}.txt", sample=SAMPLE), 
        expand("Bracken_report/{sample}.txt", sample=SAMPLE),
        expand("Contigs/flye/{sample}/assembly.fasta", sample=SAMPLE), 
        expand("Kraken_contigs/{sample}.txt", sample=SAMPLE), 
        expand("Kraken_contigs_out/{sample}.txt", sample=SAMPLE)

rule fastqc_rawreads: 
    input: 
        rawread="RawReads/{sample}.fastq.gz"
    output: 
        zip="RawQC/{sample}_fastqc.zip",
        html="RawQC/{sample}_fastqc.html"
    conda:
        "envs/fastqc_env.yaml"
    log: "logs/fastqc_{sample}.log"
    threads:
        1
    params:
        path="RawQC"
    shell:
        """
        fastqc  {input.rawread} \
                --threads {threads} \
                -o {params.path}  > {log} 2>&1 
        """ 

#rule chopper_run: 
#    input:
#        rawread="RawReads/{sample}.fastq.gz"
#    output:
#        FilteredRead="FilteredReads/{sample}.fastq.gz"
#    params:
#        threads=5
#    shell: 
#        """
#        gunzip -c {input.rawread} | \
#        chopper -q 10 -l 100 --tailcrop 10 --threads {params.threads} | \
#        gzip > {output.FilteredRead}
#        """ 

#rule fastqc_filtered: 
#    input:
#        FilteredRead="FilteredReads/{sample}.fastq.gz"
#    output:
#        zip="FilteredQC/{sample}_fastqc.zip",
#        html="FilteredQC/{sample}_fastqc.html"
#    
#    params:
#        path="FilteredQC"
#    threads:
#        5
#        
#    shell: 
#        """
#        fastqc {input.FilteredRead} --threads {threads} -o {params.path} 
 #       """

rule kraken2_run: 
    input: 
        rawread="RawReads/{sample}.fastq.gz"
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        threads=5
    output:
        output="Kraken_output/{sample}.txt",
        report="Kraken_report/{sample}.txt" 
    conda:
        "envs/kraken2_env.yaml"
    log:
        "logs/kraken2_{sample}.log"
    shell:
        """
        kraken2 --db {params.database} \
                --gzip-compressed {input.rawread} \
                --output {output.output} \
                --report {output.report} \
                --threads {params.threads} \
                --confidence 0.01 \
                --use-names   > {log} 2>&1 
        """
rule bracken_run:
    input:
        report=rules.kraken2_run.output.report
    output:
        report="Bracken_report/{sample}.txt"
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        length=100
    conda:
        "envs/bracken_env.yaml"
    log:
        "logs/bracken_{sample}.log"

    shell:
        """
        bracken -d {params.database}\
                -i {input.report}\
                -o {output.report}\
                -r {params.length}  > {log} 2>&1
        """

rule assembly_flye:
    input: 
        reads = "RawReads/{sample}.fastq.gz"
    output:
        contigs = "Contigs/flye/{sample}/assembly.fasta" 
    params:
        outdir = "Contigs/flye/{sample}",
        genome_size = "5g"
    threads:
        5
    conda:
        "envs/flye_env.yaml"
    log: 
        "logs/flye_{sample}.log"
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
        report="Kraken_contigs/{sample}.txt",
        output="Kraken_contigs_out/{sample}.txt"
    conda:
        "envs/kraken2_env.yaml"
    log:
        "logs/kraken2_Contigs_{sample}.log"
    shell:
        """
        kraken2 --db {params.database} {input.contigs} \
                --report {output.report} \
                --output {output.output} \
                --threads {params.threads} \
                --use-mpa-style \
                --use-names   > {log} 2>&1
        """

        

