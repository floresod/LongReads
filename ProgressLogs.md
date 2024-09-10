This Workflow will include multiple steps to analyze long-reads from nanopore sequencing. Hopefully with the right format.

2024-Sep-06: I am currently working on rule 1, fastqc + multiqc. Need to add multiqc code.

2024-Sep-08: I completed the fastqc and multiqc rules. I did a dry run and it seemed to work. I also started the kraken2 rule. Also made a few correctios to the kraken, bracken, flye and kraken2 on contigs rules. still need to decide which results need to be generated with snakemake. 

2024-Sep-09: Corrected the code as it was failing to run multiqc before fastqc. The problema was a bad reference to the samples when defining the wildcard SAMPLE. 


