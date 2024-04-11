Scripts used for the master thesis
# Intratumoral microorganisms and their effects on antitumor immunity

> Identificaton \
> Scripts used for the identification of intratumoral microorganisms in TCGA bam files
- 01a_combi_down_ident_loop.sh: Script reading in TCGA manifest file and submitting job for each file
- 01b_combi_down_ident.sh: Download of bam file and extraction of unaligned reads as fastq
- 01c_pipe_identification.py: Python3 script running Kraken2 and Bracken (functions in: functions_identification_pipe.py)

- 02_get_barcodes.R: R script that creates a file with UUIDs and TCGA barcodes
- 02_check_files.sh: Checks if each bam file has Bracken output

- 03_create_biom.sh: Creates BIOM file of all Bracken outputs with the 03_create_biom.py

> Analysis_OV\
> Scripts used for the different analysis, as example for ovarian cancer
- 
