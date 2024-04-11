#### !! Download of needed data !! ####

#### Input that needs to be adjusted
biom_path = "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Results_NEW/Final_nt.biom"
cancer_type = "TCGA-OV"
cancer_type_short = substr(cancer_type, 6, nchar(cancer_type))
pathout = "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"
bampath = "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_RNAseq"

## Load needed libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(phyloseq)


#### Downloads
### Downloading clinical data
clinical <- GDCquery_clinic(project = cancer_type, type = "clinical", save.csv = F)

### Downloading aliquot data
aliquot_query <- GDCquery(project = cancer_type, data.category = "Biospecimen", file.type = "xml")
# GDCdownload(aliquot_query, directory = pathout)
aliquot <- GDCprepare_clinic(aliquot_query, "aliquot", directory = pathout)

### Downloading molecular subtypes
subtypes_table <- PanCancerAtlas_subtypes()

### Downloading survival data
survival <- read.csv(paste("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2F", 
                           cancer_type_short, "_survival.txt", sep = ""),
                     sep = "\t")

# Save survival data
write.csv(survival, paste0(pathout, "/Survival_data.csv"))

### Downloading gene counts
count_query <- GDCquery(project = cancer_type, 
                        data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts")
#! GDCdownload(count_query, directory = pathout)
counts <- GDCprepare(count_query, summarizedExperiment = TRUE, directory = pathout)
count_matrix <- as.data.frame(assay(counts, "fpkm_unstrand"))
raw_matrix <- as.data.frame(assay(counts, "unstranded"))

## Saving raw and normalized count matrix
utils::write.csv(raw_matrix, file = paste0(pathout, "/counts_raw_unstrand.csv"))
utils::write.csv(count_matrix, file = paste0(pathout, "/counts_fpkm_unstrand.csv"))

#### Combine all tables by sample_id
### Loading information of needed samples
TCGA_barcodes <- read.csv(paste(bampath, "/TCGA_barcodes.txt", sep = ""), sep = " ")

patients <- sort(TCGA_barcodes$Shortened)
entity_id <- sort(TCGA_barcodes$Barcode)
samples <- sort(substr(entity_id, 1, 16))
subtype_id <- sort(substr(entity_id, 1, 20))

### Creating new data frame with ids and sample_type
new_clinical <- data.frame(submitter_id = patients,
                           sample_id = samples,
                           bcr_aliquot_barcode = entity_id,
                           subtype_id = subtype_id)

new_clinical <- new_clinical %>%
  mutate(sample_type = ifelse(substr(sample_id, 14, 15) == "01", "Primary", "Recurrent"))

## For OV (TCGA-23-1023-01R is recurrent tumor):
new_clinical <- new_clinical %>%
  mutate(sample_type = ifelse(sample_id == "TCGA-23-1023-01R", "Recurrent", sample_type))

### Merging new dataframe with clinical data
merged_dataframe <- merge(new_clinical, clinical, by = "submitter_id", all.x = T)

### Merging dataframe with aliquot based on whole entity id
merged_dataframe <- merge(merged_dataframe, aliquot, by = "bcr_aliquot_barcode")
merged_dataframe <- merged_dataframe %>% distinct(bcr_aliquot_uuid, .keep_all = TRUE)

### Merging with molecular subtypes
colnames(subtypes_table)[colnames(subtypes_table) == "pan.samplesID"] <- "subtype_id"
merged_dataframe <- merge(merged_dataframe, subtypes_table, by = "subtype_id", all.x = T)

#### Save metadata as csv
write.csv(merged_dataframe, file = paste(pathout, "/metadata.csv", sep = ""), row.names = F)

##### Loading BIOM file and adding sample data
data <- import_biom(biom_path)

## Removing _bracken_species from sample names
sample_names(data) <- substr(sample_names(data), 1, 16)

## Changing names of tax_table and removing prefixes
colnames(tax_table(data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
tax_table(data)[, colnames(tax_table(data))] <- gsub(tax_table(data)[, colnames(tax_table(data))], 
                                                     pattern = "[a-z]__", replacement = "")

## Creating sampledata object
sampledata <- sample_data(merged_dataframe)
row.names(sampledata) <- sampledata$sample_id

## Merging phyloseq object
data <- merge_phyloseq(data, sampledata)

#### Saving phyloseq objects
## With everything
saveRDS(data, file = paste(pathout, "/all.rds", sep = ""))

## Only bacteria
bacteria = subset_taxa(data, Kingdom == "Bacteria")

## Only fungi
fungi = subset_taxa(data, grepl("mycota$|Microsporidia|Nephridiophagidae", Phylum))

## Only viruses
virus = subset_taxa(data, Kingdom == "Viruses")

## All microorganisms
microorganisms = merge_phyloseq(bacteria, fungi, virus)

## Checking for samples that have under 10 reads for all three categories
empty_rows <- colSums(otu_table(microorganisms)) <=10
table(empty_rows)

## Removing those samples and saving the rds files
phyloseqs <- c("bacteria", "fungi", "virus", "microorganisms")

for (i in phyloseqs) {
  pruned_phyloseq <- prune_samples(!empty_rows, get(i))
  saveRDS(pruned_phyloseq, paste0(pathout, "/", i, ".rds"))
}
