##### !! Filtering core species after manual literature review and drawing new heatmaps !! ######
### Table of core species should have been manually reviewed by:
## 1) Checking presence in MiMeDB
## 2) Checking available literature

# Load libraries
set.seed(1234)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(pals)


##### Input #####
pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"
path_plots <- paste0(pathout, "/Plots")


##### [1] Filtering manually reviewed dataframes
## Reading in manually reviewed dataframes
bacteria <- read.csv(paste0(pathout, "/Bacteria_core.csv"))
fungi <- read.csv(paste0(pathout, "/Fungi_core.csv"))
virus <- read.csv(paste0(pathout, "/Virus_core.csv"))

## Filtering dataframes for species marked with KEEP
cores <- c("bacteria", "fungi", "virus")

for (i in cores) {
  filtered_df <- dplyr::filter(get(i), grepl("KEEP", Decision))
  df_name <- paste0(i, "_filtered")
  assign(df_name, filtered_df)
  
  tax_ids <- c()
  tax_ids <- as.character(filtered_df$Row.names)
  ids_name <- paste0(i, "_ids")
  assign(ids_name, tax_ids)
  cat("Keeping ", length(tax_ids), " species for ", i, ". \n", sep = "")
}


##### [2] Subsetting the respective phyloseqs
### Loading the phyloseqs
bacteria_phy <- readRDS(paste0(pathout, "/bacteria_core.rds"))
fungi_phy <- readRDS(paste0(pathout, "/fungi_core.rds"))
virus_phy <- readRDS(paste0(pathout, "/virus_core.rds"))

### Filtering the phyloseqs for the manually reviewed species
bacteria_fil <- prune_taxa(bacteria_ids, bacteria_phy)
saveRDS(bacteria_fil, file = paste0(pathout, "/bacteria_filtered_core.rds"))
fungi_fil <- prune_taxa(fungi_ids, fungi_phy)
saveRDS(fungi_fil, file = paste0(pathout, "/fungi_filtered_core.rds"))
virus_fil <- prune_taxa(virus_ids, virus_phy)
saveRDS(virus_fil, file = paste0(pathout, "/virus_filtered_core.rds"))

### Merging together the filtered phyloseqs
micro_fil <- merge_phyloseq(bacteria_fil, fungi_fil, virus_fil)
saveRDS(micro_fil, file = paste0(pathout, "/microorganisms_filtered_core.rds"))


###### [3] Drawing the heatmaps with manually reviewed core
### [3.1] For Bacteria
## Getting species names and ids for heatmap
bacteria_tax <- as.data.frame(tax_table(bacteria_fil))
bacteria_id <- row.names(bacteria_tax)
bacteria_species = c(paste(bacteria_tax$Genus, bacteria_tax$Species))

## Getting sample data for molecular subtypes
bacteria_sample <- as.data.frame(sample_data(bacteria_fil))

## Creating log-transformed matrix for heatmap
bacteria_log <- microbiome::transform(bacteria_fil, "log10")
bacteria_matrix <- data.matrix(as.data.frame(otu_table(bacteria_log)))

## Preparing colors for heatmap
col_heatmap <- rev(ocean.deep(n = 255))

## Row annotation
# Getting colors for each phyla 
bacteria_unique_phyla <- unique(bacteria_tax$Phylum)
bacteria_num_phyla <- length(bacteria_unique_phyla)

bacteria_colors_rows <- cubicl(n = bacteria_num_phyla)
bacteria_color_mapping_rows <- setNames(bacteria_colors_rows, bacteria_unique_phyla)

# Creating row annotation
bacteria_row_ha <- rowAnnotation(Phylum = bacteria_tax$Phylum, 
                        col = list(Phylum = bacteria_color_mapping_rows),
                        gp = gpar(col = "black"),
                        simple_anno_size = unit(0.3, "cm"),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

## Column annotation: molecular subtype 
# Getting colors for each molecular subtype
bacteria_unique_subtypes <- unique(bacteria_sample$Subtype_mRNA)
bacteria_num_subtypes <- length(bacteria_unique_subtypes)

bacteria_colors_cols <- trubetskoy(n = bacteria_num_subtypes)
bacteria_color_mapping_cols <- setNames(bacteria_colors_cols, bacteria_unique_subtypes)

# Removing NA from named vector
bacteria_color_mapping_cols <- bacteria_color_mapping_cols[!is.na(names(bacteria_color_mapping_cols))]

# Generating column annotation
bacteria_col_ha <- columnAnnotation(Subtype = bacteria_sample$Subtype_mRNA, 
                        col = list(Subtype = bacteria_color_mapping_cols),
                        simple_anno_size = unit(0.3, "cm"),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

# Saving heatmap for bacteria
png(file=paste0(path_plots, "/Bacteria_Heatmap_filtered.png"),
    width=(nsamples(bacteria_fil) * 0.08), height=(ntaxa(bacteria_fil) * 0.4), units="cm",res=1200)
ht_bac <- Heatmap(bacteria_matrix, 
              show_column_names = FALSE, row_labels = bacteria_species,
              row_names_gp = gpar(fontsize = 9, fontface = "italic"),
              col = col_heatmap,
              name = "log10",
              column_title = "Samples",
              column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 10),
              right_annotation = bacteria_row_ha,
              top_annotation = bacteria_col_ha)
draw(ht_bac, merge_legend = TRUE)
dev.off()

### [3.2] For fungi
## Getting species names and ids for heatmap
fungi_tax <- as.data.frame(tax_table(fungi_fil))
fungi_id <- row.names(fungi_tax)
fungi_species = c(paste(fungi_tax$Genus, fungi_tax$Species))

## Getting sample data for molecular subtypes
fungi_sample <- as.data.frame(sample_data(fungi_fil))

## Creating log-transformed matrix for heatmap
fungi_log <- microbiome::transform(fungi_fil, "log10")
fungi_matrix <- data.matrix(as.data.frame(otu_table(fungi_log)))

## Preparing colors for heatmap
col_heatmap <- rev(ocean.deep(n = 255))

## Row annotation
# Getting colors for each phyla 
fungi_unique_phyla <- unique(fungi_tax$Phylum)
fungi_num_phyla <- length(fungi_unique_phyla)

fungi_colors_rows <- cubicl(n = fungi_num_phyla)
fungi_color_mapping_rows <- setNames(fungi_colors_rows, fungi_unique_phyla)

# Creating row annotation
fungi_row_ha <- rowAnnotation(Phylum = fungi_tax$Phylum, 
                                 col = list(Phylum = fungi_color_mapping_rows),
                                 gp = gpar(col = "black"),
                                 simple_anno_size = unit(0.3, "cm"),
                                 annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

## Column annotation: molecular subtype 
# Getting colors for each molecular subtype
fungi_unique_subtypes <- unique(fungi_sample$Subtype_mRNA)
fungi_num_subtypes <- length(fungi_unique_subtypes)

fungi_colors_cols <- trubetskoy(n = fungi_num_subtypes)
fungi_color_mapping_cols <- setNames(fungi_colors_cols, fungi_unique_subtypes)

# Removing NA from named vector
fungi_color_mapping_cols <- fungi_color_mapping_cols[!is.na(names(fungi_color_mapping_cols))]

# Generating column annotation
fungi_col_ha <- columnAnnotation(Subtype = fungi_sample$Subtype_mRNA, 
                                    col = list(Subtype = fungi_color_mapping_cols),
                                    simple_anno_size = unit(0.3, "cm"),
                                    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

# Saving heatmap for fungi
png(file=paste0(path_plots, "/Fungi_Heatmap_filtered.png"),
    width=(nsamples(fungi_fil) * 0.08), height=(ntaxa(fungi_fil) * 0.4), units="cm",res=1200)
ht_fun <- Heatmap(fungi_matrix, 
                  show_column_names = FALSE, row_labels = fungi_species,
                  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
                  col = col_heatmap,
                  name = "log10",
                  column_title = "Samples",
                  column_title_side = "bottom",
                  column_title_gp = gpar(fontsize = 10),
                  right_annotation = fungi_row_ha,
                  top_annotation = fungi_col_ha)
draw(ht_fun, merge_legend = TRUE)
dev.off()

### [3.3] For viruses
## Getting species names and ids for heatmap
virus_tax <- as.data.frame(tax_table(virus_fil))
virus_id <- row.names(virus_tax)
virus_species = c(paste(virus_tax$Genus, virus_tax$Species))

## Getting sample data for molecular subtypes
virus_sample <- as.data.frame(sample_data(virus_fil))

## Creating log-transformed matrix for heatmap
virus_log <- microbiome::transform(virus_fil, "log10")
virus_matrix <- data.matrix(as.data.frame(otu_table(virus_log)))

## Preparing colors for heatmap
col_heatmap <- rev(ocean.deep(n = 255))

## Row annotation
# Getting colors for each phyla 
virus_unique_phyla <- unique(virus_tax$Phylum)
virus_num_phyla <- length(virus_unique_phyla)

virus_colors_rows <- cubicl(n = virus_num_phyla)
virus_color_mapping_rows <- setNames(virus_colors_rows, virus_unique_phyla)

# Creating row annotation
virus_row_ha <- rowAnnotation(Phylum = virus_tax$Phylum, 
                              col = list(Phylum = virus_color_mapping_rows),
                              gp = gpar(col = "black"),
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

## Column annotation: molecular subtype 
# Getting colors for each molecular subtype
virus_unique_subtypes <- unique(virus_sample$Subtype_mRNA)
virus_num_subtypes <- length(virus_unique_subtypes)

virus_colors_cols <- trubetskoy(n = virus_num_subtypes)
virus_color_mapping_cols <- setNames(virus_colors_cols, virus_unique_subtypes)

# Removing NA from named vector
virus_color_mapping_cols <- virus_color_mapping_cols[!is.na(names(virus_color_mapping_cols))]

# Generating column annotation
virus_col_ha <- columnAnnotation(Subtype = virus_sample$Subtype_mRNA, 
                                 col = list(Subtype = virus_color_mapping_cols),
                                 simple_anno_size = unit(0.3, "cm"),
                                 annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

# Saving heatmap for virus
png(file=paste0(path_plots, "/Virus_Heatmap_filtered.png"),
    width=(nsamples(virus_fil) * 0.08), height=(ntaxa(virus_fil) * 0.4), units="cm",res=1200)
ht_vir <- Heatmap(virus_matrix, 
                  show_column_names = FALSE, row_labels = virus_species,
                  row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                  col = col_heatmap,
                  name = "log10",
                  column_title = "Samples",
                  column_title_side = "bottom",
                  column_title_gp = gpar(fontsize = 10),
                  right_annotation = virus_row_ha,
                  top_annotation = virus_col_ha)
draw(ht_vir, merge_legend = TRUE)
dev.off()






