This Workflow will include multiple steps to analyze long-reads from nanopore sequencing. Hopefully with the right format.

2024-Sep-06: I am currently working on rule 1, fastqc + multiqc. Need to add multiqc code.

2024-Sep-08: I completed the fastqc and multiqc rules. I did a dry run and it seemed to work. I also started the kraken2 rule. Also made a few correctios to the kraken, bracken, flye and kraken2 on contigs rules. still need to decide which results need to be generated with snakemake. 

2024-Sep-09: Corrected the code as it was failing to run multiqc before fastqc. The problema was a bad reference to the samples when defining the wildcard SAMPLE. 

2024-Sep-10: Found that the problem was that the fastqc.gz files were not compressed, they were actually fastqc only. Changing the name fixed the problem. There is a problem with multiqc. I will leave it out for now. **Need to correct this in the future**

2024-Sep-11: Added Kraken2, Bracken, and Flye code to the snakemake file. It works well and produces the contigs. 

2024-Sep-12: Added the Medaka code. It works and generates the consensus assemblies.

2024-Oct-04: Added code to move contigs generated with medaka to results/FinalContigs. 

2024-Oct-07: Added and tested checkm2 code. 

2024-Oct-10: Added gtdbtk code, no tested yet
