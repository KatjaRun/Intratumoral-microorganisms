Scripts from the master thesis:
# Intratumoral microorganisms and their effects on antitumor immunity

> Identificaton \
> Identification of intratumoral microorganisms in TCGA bam files
- 01a_combi_down_ident_loop.sh: Script reading in TCGA manifest file and submitting job for each file
- 01b_combi_down_ident.sh: Download of bam file and extraction of unaligned reads as fastq
- 01c_pipe_identification.py: Python3 script running Kraken2 and Bracken (functions in: functions_identification_pipe.py)

- 02_get_barcodes.R: R script that creates a file with UUIDs and TCGA barcodes
- 02_check_files.sh: Checks if each bam file has Bracken output

- 03_create_biom.sh: Creates BIOM file of all Bracken outputs with the 03_create_biom.py


> Analysis_OV\
> Different analysis performed on each cancer type, as example for ovarian cancer
- 01_combine_metadata_biom.R: Downloading needed metadata and adds it to the BIOM file as phyloseq object
- 02_Bacteria/Fungi/Viruses.R: Decontamination and identification of core species for each kingdom, respectively
  
3. Manual literature review of core species

- 04_Filtering_manual_review.R: Filtering of species based on manual literature review
- 05_SPIEC_EASI.R: Assessing microbial communities through co-occurrence
- 06_Correlation_CIBER.R: Calculating immune cell fractions and then performing correlation with core species abundance
- 07_Survival_Analysis.R: Performing Cox proportional hazards model to assess impact on overall survival of core species
- 08_Correaltion_GSVA.R: Computing GSVA scores and correlating them to core species abundance
- 09_Interesting_species.R: Creating an Euler diagram to assess species significant in performed analyses

> Additional
- DESeq_PAAD_NAT.R: Performing differential abundance analysis for pancreatic tumor and normal samples
- PCoA_Beta.R: Calculating Bray distance for ovarian and pancreatic adenocarcinoma

- extract_int_reads.R: Extracts reads from Kraken2 files based on TaxID
- STAR_int_reads.sh: Aligns extracted reads to coding sequences for the respective species (TaxID)
- Depth_Plots.R: Calculating RPKM values per coding sequence
