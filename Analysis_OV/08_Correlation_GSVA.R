#### !! Linear regression with GSVA immune parameters !! #####

pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"

# Libraries
library(GSVA)
library(biomaRt)
library(phyloseq)
library(tidyverse)
library(venneuler)

###### [1] Performing GSVA for several immune parameters ######
## Expression data table
counts_fpkm <- as.data.frame(read.csv(paste0(pathout, "/counts_fpkm_unstrand.csv")))

## Stripping the version number of Ensembl Gene IDs
counts_fpkm$X <- sub("\\..*", "", counts_fpkm$X)

## Getting gene symbol as row.names and removing Ensembl Gene IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)

## Merging the dataframes
counts <- merge(conversion, counts_fpkm, by.x = "ensembl_gene_id", by.y = "X")

## Putting gene names as row.names and removing the excessive columns
row.names(counts) <- make.unique(counts$external_gene_name)
data <- subset(counts, select = -c(ensembl_gene_id, external_gene_name))

## Getting same sample names as in the microbiome data
colnames(data) <- gsub("\\.", "-", colnames(data))
colnames(data) <- substr(colnames(data), 1, 16)

## Preparing gene sets for GSVA
gs <- read.csv(paste0(pathout, "/GSVA_gene_sets/Gene_sets.csv"), header = T)
gs_list = split(gs$Gene, gs$Annotation)

### Performing GSVA 
gsva_res <- as.data.frame(GSVA::gsva(as.matrix(data), gs_list, method = "gsva", kcdf = "Gaussian"))
write.csv(gsva_res, paste0(pathout, "/GSVA_results.csv"))

##### [2] Reading in abundance data #####
core <- readRDS(paste0(pathout, "/microorganisms_filtered_core.rds"))
core_clr <- microbiome::transform(core, "clr")

# Getting tax table with species names
core_tax <- as.data.frame(tax_table(core_clr))
core_tax$Names <-  ifelse(core_tax$Kingdom %in% c("Bacteria", "Eukaryota"), 
                          paste(core_tax$Genus, core_tax$Species, sep = " "),
                          ifelse(core_tax$Kingdom == "Viruses",
                                 ifelse(nchar(core_tax$Species) == 1, 
                                        paste(core_tax$Genus, core_tax$Species, sep = " "), 
                                        core_tax$Species),
                                 "Unknown_Kingdom"))

# Saving OTU table and changing the rownames to species name
core_OTU <- as.data.frame(otu_table(core_clr))
row.names(core_OTU) <- core_tax$Names

## Filtering out samples not in the OTU table
gsva_res <- gsva_res[ , colnames(gsva_res) %in% colnames(core_OTU)]

### Combining OTU_table and GSVA data
data_lm <- rbind(gsva_res, core_OTU)

data_lm <- as.data.frame(t(data_lm))

### Filtering out outliers by three times standard deviation
columns <- colnames(data_lm)

for (i in columns) {
  mean <- mean(data_lm[,i])
  stbw <- sd(data_lm[,i])
  threshold <- 2
  
  clean_data_lm <- data_lm[!(data_lm[,i] < (mean - threshold * stbw) | 
                         data_lm[,i] > (mean + threshold * stbw)), ]
}


###### [3] Calculating regression ######
# Function to get p-value (https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Looping parameters
species <- row.names(core_OTU)
gene_sets <- unique(gs$Annotation)

## Preparing file for saving
header <- paste("Gene_set", "Species", "LM_R2", "LM_p", "Corr_r", "Corr_p", sep = "\t")
write.table(header, file= file.path(pathout, "/GSVA_gene_sets/Regression_table.txt"), append = TRUE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

## Looping through all parameters and writing to a file
for (i in gene_sets) {
  for (j in species) {
    linear_model <- lm(clean_data_lm[,j] ~ clean_data_lm[,i], data = clean_data_lm) # Regression
    r_squared <- summary(linear_model)$r.squared # Extracting R²
    lm_p <- lmp(linear_model) # Extracting p
    
    # Performing correlation analysis
    correl <- rcorr(clean_data_lm[,i], clean_data_lm[,j], "pearson")
    corr_r <- correl$r[1, 2]
    corr_p <- correl$P[1, 2]
    
    # Appending result to table
    text <-paste(i, j, r_squared, lm_p, corr_r, corr_p, sep="\t")
    write.table(text, file= file.path(pathout, "/GSVA_gene_sets/Regression_table.txt"),append = TRUE, quote = FALSE, sep = "\t", 
                eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  }
}

## Reading in table
regression <- read.table(file.path(pathout, "/GSVA_gene_sets/Regression_table.txt"), header = TRUE, sep = "\t")

## FDR correction
regression$Corr_adj_p <- p.adjust(regression$Corr_p, "fdr")

# Filtering for significant immune cells
regression_fil <- regression[regression$Corr_adj_p < 0.05, ]

write.csv(regression_fil, file = paste0(pathout, "/GSVA_gene_sets/Regression_table_fil.csv"))

##### [5] Drawing heatmap for correlations #####
## Performing correlation without outliers
corr_matrix <- rcorr(as.matrix(clean_data_lm), type = "pearson")

## Filtering the dataframe to get species in rows and immune parameters in columns
# Species with minimal one significant parameter
sig_species <- unique(regression_fil$Species)
corr_r <- as.data.frame(corr_matrix$r)
corr_ht <- corr_r[row.names(corr_r) %in% sig_species, colnames(corr_r) %in% gene_sets]

# Extracting species order
sig_species_order <- row.names(corr_ht)

## Getting p values and filtering for same values
corr_p <- as.data.frame(corr_matrix$P)
corr_p_fil <- corr_p[row.names(corr_p) %in% sig_species, colnames(corr_p) %in% gene_sets]

## Adjusting p values
corr_adjp <- matrix(p.adjust(as.vector(as.matrix(corr_p_fil)), method='fdr'),ncol=ncol(corr_p_fil))


### Drawing heatmap for species with immune cells
## Creating row annotations
# Tax table 
sig_tax <- core_tax[core_tax$Names %in% sig_species_order, ]
sig_tax$Names <- factor(sig_tax$Names, levels = sig_species_order)
sig_tax <- sig_tax[order(sig_tax$Names), ]

# Replacing Eukaryota with fungi
sig_tax$Kingdom <- ifelse(sig_tax$Kingdom == "Eukaryota", "Fungi", sig_tax$Kingdom)


# Saving tax table of significant species
write.csv(sig_tax, paste0(pathout, "/GSVA_species.csv"))

# Replacing odd species names
rownames(corr_ht)[rownames(corr_ht) == "haeme3"] <- "Gemykibivirus haeme3"
rownames(corr_ht)[rownames(corr_ht) == "Pseudomonas uncultured Pseudomonas sp."] <- "Pseudomonas sp."
rownames(corr_ht)[rownames(corr_ht) == "Desulfosudis Deltaproteobacteria bacterium"] <- "Desulfosudis sp."
rownames(corr_ht)[rownames(corr_ht) == "Ambispora uncultured Glomeromycotina"] <- "Glomeromycotina sp."
rownames(corr_ht)[rownames(corr_ht) == "Ambispora uncultured Glomeromycotina"] <- "Glomeromycotina sp."

# And changing column names
colnames(corr_ht)[colnames(corr_ht) == "IFNG_Ayers"] <- "IFNG (Ayers)"
colnames(corr_ht)[colnames(corr_ht) == "IL10_pathway"] <- "IL-10 pathway"
colnames(corr_ht)[colnames(corr_ht) == "IL12_pathway"] <- "IL-12 pathway"
colnames(corr_ht)[colnames(corr_ht) == "IL4_pathway"] <- "IL-4 pathway"
colnames(corr_ht)[colnames(corr_ht) == "IL6_pathway"] <- "IL-6 pathway"
colnames(corr_ht)[colnames(corr_ht) == "KEGG_CGAS_STING"] <- "CGAS-STING pathway"
colnames(corr_ht)[colnames(corr_ht) == "KEGG_JAK_STAT"] <- "JAK-STAT"
colnames(corr_ht)[colnames(corr_ht) == "NFKB_pathway"] <- "NFKB pathway"
colnames(corr_ht)[colnames(corr_ht) == "PID_IL8_CXCR2_pathway"] <- "IL-8 CXCR2"
colnames(corr_ht)[colnames(corr_ht) == "RAS_pathway"] <- "RAS pathway"
colnames(corr_ht)[colnames(corr_ht) == "TGFB_pathway"] <- "TGFB pathway"
colnames(corr_ht)[colnames(corr_ht) == "TNFR1_pathway"] <- "TNFR1 pathway"
colnames(corr_ht)[colnames(corr_ht) == "VEGF_pathway"] <- "VEGF pathway"
colnames(corr_ht)[colnames(corr_ht) == "GOBP_APP"] <- "Antigen pres. and proc."
colnames(corr_ht)[colnames(corr_ht) == "CYT"] <- "Cytolytic Activity"

# For kingdoms
color_mapping_kingdoms <- setNames(c("#377eb8", "#984ea3", "#e41a1c"), c("Bacteria", "Fungi", "Viruses"))

col_Kingdom <- rowAnnotation(Kingdom = sig_tax$Kingdom, 
                             col = list(Kingdom = color_mapping_kingdoms),
                             gp = gpar(col = "black"),
                             simple_anno_size = unit(0.3, "cm"),
                             annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                             annotation_legend_param = list(labels_gp = gpar(fontsize = 17),
                                                            title_gp = gpar(fontsize = 17, fontface = "bold")))

# For phyla
unique_phyla <- unique(sig_tax$Phylum)
num_phyla <- length(unique_phyla)

colors_phyla <- trubetskoy(n = num_phyla)
color_mapping_phyla <- setNames(colors_phyla, unique_phyla)

col_Phyla <- rowAnnotation(Phylum = sig_tax$Phylum, 
                           col = list(Phylum = color_mapping_phyla),
                           gp = gpar(col = "black"),
                           simple_anno_size = unit(0.3, "cm"),
                           annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                           annotation_legend_param = list(labels_gp = gpar(fontsize = 17),
                                                          title_gp = gpar(fontsize = 17, fontface = "bold")))


# Drawing heatmap
ht <- Heatmap(as.matrix(corr_ht),
              name = "Pearson \ncorrelation",
              col = colorRamp2(c(-0.5, -0.25, 0, 0.25, 0.5), c("#67001f", "#b2182b", "white", "#2166ac", "#053061")),
              right_annotation = c(col_Kingdom, col_Phyla),
              row_km = 2,
              row_names_gp = gpar(fontsize = 16, fontface = "italic"),
              column_names_gp = gpar(fontsize = 16),
              heatmap_legend_param = list(labels_gp = gpar(fontsize = 17),
                                          title_gp = gpar(fontsize = 17, fontface = "bold")),
              cell_fun = function(j, i, x, y, w, h, fill) {
                asterisk_position_y <- y - h/2.5
                
                ifelse(corr_adjp[i, j] < 0.001, grid.text("***", x, asterisk_position_y, gp = gpar(fontsize = 15)),
                       ifelse(corr_adjp[i, j] < 0.01, grid.text("**", x, asterisk_position_y, gp = gpar(fontsize = 15)),
                              ifelse(corr_adjp[i, j] < 0.05, grid.text("*", x, asterisk_position_y, gp = gpar(fontsize = 15)), 
                                     grid.text(""))))
              }
)

# Saving
png(paste0(pathout, "/Plots/GSVA_species_heatmap.png"), width=30,height=35,units="cm",res=1200)
ht_draw <- draw(ht, heatmap_legend_side = "left", merge_legend = T, padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

### Heatmap for only immune parameters
corr_immune <- corr_r[row.names(corr_r) %in% gene_sets, colnames(corr_r) %in% gene_sets]

# Adjusting p values
corr_immune_p <- corr_p[row.names(corr_p) %in% gene_sets, colnames(corr_p) %in% gene_sets]
corr_immune_padj <- matrix(p.adjust(as.vector(as.matrix(corr_immune_p)), method='bonferroni'),ncol=ncol(corr_immune_p))

# Correlation values
ht_immune <- Heatmap(as.matrix(corr_immune),
              name = "Pearson \ncorrelation",
              col = colorRamp2(c(-0.5, 0, 0.5, 1), c("#67001f", "white", "#2166ac", "#053061")),
              cell_fun = function(j, i, x, y, w, h, fill) {
                asterisk_position_y <- y - h/2.5
                
                ifelse(corr_immune_padj[i, j] < 0.001, grid.text("***", x, asterisk_position_y, gp = gpar(fontsize = 15)),
                       ifelse(corr_immune_padj[i, j] < 0.01, grid.text("**", x, asterisk_position_y, gp = gpar(fontsize = 15)),
                              ifelse(corr_immune_padj[i, j] < 0.05, grid.text("*", x, asterisk_position_y, gp = gpar(fontsize = 15)), 
                                     grid.text(""))))
              }
)

# Saving
png(paste0(pathout, "/Plots/GSVA_GSVA_heatmap.png"), width=20,height=20,units="cm",res=1200)
ht_draw_immune <- draw(ht_immune, heatmap_legend_side = "left", merge_legend = T, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()



