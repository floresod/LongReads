##################################
###### LONG READS ANALYSES #######
##################################


################################
#### Loado python libraries ####
################################
import glob
#import os 


##########################
#### Global Variables ####
##########################
SAMPLE,=glob_wildcards("../resources/Data/Fastq/{sample}.fastq")

########################
#### Global outputs ####
########################
rule all:
    input:
        expand("../resources/Outputs/fastqc_rawreads/{sample}_fastqc.html", sample=SAMPLE), 
        expand("../results/kraken2_rr/{sample}.txt", sample=SAMPLE),
        #expand("../results/bracken/{sample}.txt", sample=SAMPLE),
        expand("../resources/Outputs/flye/{sample}/assembly.fasta", sample=SAMPLE),
        expand("../resources/Outputs/medaka/{sample}/consensus.fasta", sample=SAMPLE), 
        expand("../results/FinalContigs/{sample}.fasta", sample = SAMPLE),
        #expand("../results/checkm2/{sample}.tsv", sample=SAMPLE),
       # expand("../resources/Outputs/gtdbtk/{sample}/", sample=SAMPLE),
        "../results/checkm2/checkm2_report.tsv",
        #expand("../results/card/{sample}.tsv", sample=SAMPLE),
       # expand("../results/mobileOG/{sample}.tsv", sample=SAMPLE),
       # expand("../results/kraken2_contigs/{sample}_output.tsv", sample=SAMPLE),
        expand("../results/kraken2_contigs/{sample}.tsv", sample=SAMPLE)

###################################
#### QQuality Control Raw Reads ####
###################################
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

#########################################
#### Taxonomic Classificationi Reads ####
#########################################
rule kraken2_rr:
    input:
        rawread="../resources/Data/Fastq/{sample}.fastq"
    output:
        report="../results/kraken2_rr/{sample}.txt", 
        output="../resources/Outputs/kraken2_rr/{sample}.txt"
    params:
        database="../../../Databases/k2_standard_08gb_20240904", 
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

##################
#### Assembly ####
##################
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
                         --bacteria \
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

#########################
#### Quality Contigs ####
#########################
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
        
rule combine_cm2reports:
    input:
        tsv = expand("../results/checkm2/{sample}.tsv", sample=SAMPLE)
    output:
        comb_report="../results/checkm2/checkm2_report.tsv"
    log:
        "../resources/Logs/combined_cmreport/log.log"
    shell:
        """
         # Write the header from the first file to the combined file
        head -n 1 {input.tsv[0]} > {output.comb_report}

        # Append the content of each TSV, skipping the header in subsequent files
        for tsv in {input.tsv}; do
            tail -n +2 $tsv >> {output.comb_report}
       done
        """

##########################################
#### Taxonomic Classification Contigs ####
##########################################
rule gtdbtk_classify:
    input:
        contigs = "../resources/Outputs/medaka/{sample}"
    output:
        report_dir = directory("../resources/Outputs/gtdbtk/{sample}")
        #combined_report = "../results/ClassifiedContigs/gtdbtk_all_reports.tsv"
    log:
        "../resources/Logs/gtdbtk/{sample}.log"
    resources:
        threads = 30
    conda:
        "../envs/gtdbtk_env.yaml"
    shell:
        """
        gtdbtk classify_wf  --genome_dir {input.contigs} \
                            --out_dir {output.report_dir} \
                            --skip_ani_screen \
                            -x fasta \
                            --cpus {resources.threads} \
                            --pplacer_cpus {resources.threads} > {log} 2>&1
        """

#######################
#### CARD Aligment ####
#######################
rule diamond_card:
    input:
        contigs="../resources/Outputs/medaka/{sample}/consensus.fasta",
        database="../../../Databases/card/card_v3.3.0.dmnd"
    output:
        "../results/card/{sample}.tsv"
    log:
        "../resources/Logs/card/{sample}.log"
    conda:
        "../envs/diamond_env.yaml"
    shell:
        """
        diamond blastx -d {input.database} -q {input.contigs} \
                -o {output} --id 95 --subject-cover 90 \
                --ultra-sensitive > {log} 2>&1
        """

##########################
#### MobileOG Aligment ###
##########################
rule diamond_mobileOG:
    input:
        contigs="../resources/Outputs/medaka/{sample}/consensus.fasta",
        database="../../../Databases/mobileOG/mobileOG.dmnd"
    output:
        "../results/mobileOG/{sample}.tsv"
    log:
        "../resources/Logs/mobileOG/{sample}.log"
    conda:
        "../envs/diamond_env.yaml"
    shell:
        """
        diamond blastx -d {input.database} -q {input.contigs} \
                    -o {output} --id 95 --subject-cover 95 \
                    --ultra-sensitive > {log} 2>&1
        """

#######################################################
#### Taxonomic Classification of Contigs - Kraken2 ####
#######################################################
rule kraken2_contigs:
    input:
        contigs="../resources/Outputs/medaka/{sample}/consensus.fasta"
    output:
        report="../resources/Outputs/kraken2_contigs/{sample}_report", 
        output="../resources/Outputs/kraken2_contigs/{sample}_output.tsv"
    params:
        database="../../../Databases/k2_standard_08gb_20240605", 
        threads="6"
    log:
        "../resources/Logs/kraken2_contigs/{sample}.log"
    conda: 
        "../envs/kraken2_env.yaml"
    shell:
        """
        kraken2 --db {params.database} \
                {input.contigs} \
                --report {output.report} \
                --report-minimizer-data \
                --output {output.output} \
                --threads {threads} \
                --confidence 0.0005 \
                --use-names > {log} 2>&1
        """
rule contigs_taxa:
    input:
        kk2_output="../resources/Outputs/kraken2_contigs/{sample}_output.tsv"
    output:
        classified_contigs="../results/kraken2_contigs/{sample}.tsv"
    log:
        "../resources/Logs/kraken2_contig_outputs/{sample}.log"
    shell:
        """
        #create directory
        mkdir -p ../results/kraken2_contigs

        cut -f 2-3 {input.kk2_output} > {output.classified_contigs}
        """
        
        

 


