##### !! Pipeline for microorganisms !! ####

#### Libraries
set.seed(43)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(MASS)
library(pals)

#### Input
pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"
path_plots <- paste0(pathout, "/Plots")

##### [1] Loading and putting together decontaminated phyloseq objects
bacteria_decontam <- readRDS(paste0(pathout, "/bacteria_decontam.rds"))
fungi_decontam <- readRDS(paste0(pathout, "/fungi_decontam.rds"))
virus_decontam <- readRDS(paste0(pathout, "/virus_decontam.rds"))

## Merging every phyloseq
micro_decontam <- merge_phyloseq(bacteria_decontam, fungi_decontam, virus_decontam)

## Saving decontaminated microorganisms
saveRDS(micro_decontam, paste0(pathout, "/microorganisms_decontam.rds"))

##### [2] Analysis of alpha diversity

### Creating table with diversity parameters
## Should be used with the count data
# Calculating diversity, richness and rarity of the samples
diversity <- microbiome::diversity(micro_decontam)
richness <- microbiome::richness(micro_decontam)
rarity <- microbiome::rarity(micro_decontam)
dominance <- microbiome::dominance(micro_decontam)

# Combining all data frames and adding a constant x-axis for plotting
micro_alpha <- cbind(diversity, richness, rarity, dominance)
micro_alpha$constant_x <- 0

## Creating ViolinPLot for observed species richness
observed <- ggplot(micro_alpha, aes(x = constant_x, y = observed)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Observed species richness") +
  ylim(0, 10000) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for species richness (Chao1)
chao1 <- ggplot(micro_alpha, aes(x = constant_x, y = chao1)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("chao1 species richness") +
  ylim(0, 10000) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for Shannon diversity
shannon <- ggplot(micro_alpha, aes(x = constant_x, y = shannon)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Shannon diversity") +
  ylim(0, 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for inverse simpson
inv_simpson <- ggplot(micro_alpha, aes(x = constant_x, y = inverse_simpson)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Inverse Simpson Index") +
  ylim(0, 100) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for dominance core abundance
dominance_core <- ggplot(micro_alpha, aes(x = constant_x, y = core_abundance * 100)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Dominance of core abundance [%]") +
  ylim(0, 100) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for low abundance
low_abundance <- ggplot(micro_alpha, aes(x = constant_x, y = low_abundance * 100)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Species with low abundance [%]") +
  ylim(0, 100) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Combining six plots
ggarrange(observed, chao1,
          shannon, inv_simpson, 
          dominance_core, low_abundance,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)

ggsave(paste(path_plots, "/Micro_alpha.png", sep = ""),
       width = 22, height = 25.5, units = "cm")


##### [3] Analysis of beta diversity
### Ordination using NMDS and bray distance
ord_NMDs <- phyloseq::ordinate(micro_decontam, "NMDS", "bray", maxit=1000)

options(repr.plot.width = 3, repr.plot.height = 3)
plot_ordination(micro_decontam, ord_NMDs, type = "samples", color = "Subtype_mRNA", shape = "sample_type",
                title = "Microbes: NMDs of samples") +
  theme_light() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1)

ggsave(paste0(path_plots, "/Micro_NMDS.png"), height = 10, width = 15, unit = "cm")

### Ordination using PCoA and bray distance
ord_PCOA <- ordinate(micro_decontam, "PCoA")

options(repr.plot.width = 3, repr.plot.height = 3)
plot_ordination(micro_decontam, ord_PCOA, type = "samples", color = "Subtype_mRNA", shape = "sample_type",
                title = "Microbes: PCoA of samples") +
  theme_light() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1)

ggsave(paste0(path_plots, "/Micro_PCoA.png"), height = 10, width = 15, unit = "cm")


##### [4] Multilabel feature selection

### Preparing the data
## Filtering out samples with NA as subtype
boo <- c(is.na(as.data.frame(sample_data(fungi_decontam)$Subtype_mRNA)))
micro_subtypes <- prune_samples(!boo , fungi_decontam)

## Getting subtypes as variables
subtypes <- c("Mesenchymal", "Proliferative", "Differential", "Immunoreactive")
subtypes_data <- as.factor(sample_data(micro_subtypes)$Subtype_mRNA)

## Getting OTU table as feature data 
feature_data <- as.data.frame(otu_table(micro_subtypes))
feature_data <- as.data.frame(t(feature_data))
feature_data$Subtype_mRNA <- as.factor(sample_data(micro_subtypes)$Subtype_mRNA)

## Creating 4 new columns which represent the subtypes
feature_data <- feature_data %>%
  mutate(
    Mesenchymal = if_else(Subtype_mRNA == "Mesenchymal", 1, 0),
    Proliferative = if_else(Subtype_mRNA == "Proliferative", 1, 0),
    Differential = if_else(Subtype_mRNA == "Differential", 1, 0),
    Immunoreactive = if_else(Subtype_mRNA == "Immunoreactive", 1, 0)
  )

### Looping through each variable and performing logistic regression
top_features <- c()

for (i in subtypes) {
  formel <- as.formula(paste(i, "~ ."))
  logistic_model <- glm(formel, data = feature_data, family = binomial)
  print(summary(logistic_model))
  current_df <- as.data.frame(coef(logistic_model))
  current_df$tax_id <- gsub("`", "", row.names(current_df))
  current_df <- current_df[!(row.names(current_df) == "(Intercept)"), ]
  current_df <- current_df[order(-coef(logistic_model)), ]
  assign(paste0(i, "_features"), current_df)
  top_features <- c(top_features, as.character(current_df$tax_id[1:100]))
}

# Subset the original phyloseq object to keep only the selected features
fungi_subset <- prune_taxa(top_features, fungi_decontam)

### Ordination using PCoA and bray distance
ord_PCOA <- ordinate(fungi_subset, "PCoA", "bray")

options(repr.plot.width = 3, repr.plot.height = 3)
plot_ordination(fungi_subset, ord_PCOA, type = "samples", color = "Subtype_mRNA", shape = "sample_type",
                title = "Fungi (MFS): PCoA of samples") +
  theme_light() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1)

ggsave(paste0(path_plots, "/Fungi_MFS_PCoA.png"), height = 10, width = 15, unit = "cm")


##### [5] Microorganisms core

### Defining core microbial species (>0.0001 relative abundance and present in 80% of all samples)
## Tax glom at species level
bacteria_spe <- readRDS(paste0(pathout, "/bacteria_species.rds"))
fungi_spe <- readRDS(paste0(pathout, "/fungi_species.rds"))
virus_spe <- readRDS(paste0(pathout, "/virus_species.rds"))

micro_spe = merge_phyloseq(bacteria_spe, fungi_spe, virus_spe)

saveRDS(micro_spe, paste0(pathout, "/microorganisms_species.rds"))

## Transforming data into relative abundance and setting parameters for filtering
micro_spe_rel = microbiome::transform(micro_spe, "compositional")
prev = 0.8
abund = (1/ntaxa(micro_spe) * 2)

# Getting mean species abundance
total_abundance <- rowSums(otu_table(micro_spe_rel), na.rm = TRUE)
cat("Mean total abundance: ", mean(total_abundance, na.rm = TRUE))
cat("Median total abundance: ", median(total_abundance, na.rm = TRUE))

## Defining the Core species
micro_core = core(micro_spe_rel, detection = abund, prevalence = prev, include.lowest = T)

## Getting species names and ids for heatmap
taxonomy_core <- as.data.frame(tax_table(micro_core))
core_id <- row.names(taxonomy_core)
core_species = c(paste(taxonomy_core$Genus, taxonomy_core$Species))

## Getting sample data for molecular subtypes
sample_core <- as.data.frame(sample_data(micro_core))

## Preparing matrix for heatmap
micro_core_prune = prune_taxa(core_id, micro_spe)

# Saving untransformed data of core species
saveRDS(micro_core_prune, paste0(pathout, "/micro_core.rds"))

## Creating log-transformed matrix for heatmap
micro_core_log = microbiome::transform(micro_core_prune, "log10")
core_matrix <- data.matrix(as.data.frame(otu_table(micro_core_log)))

## Preparing colors for heatmap
col_heatmap <- rev(cividis(n = 255))

## Row annotation
# Getting colors for each Kingdom
unique_phyla <- unique(taxonomy_core$Kingdom)
num_phyla <- length(unique_phyla)

colors_rows <- cubicl(n = num_phyla)
color_mapping_rows <- setNames(colors_rows, unique_phyla)

# Creating row annotation
row_ha <- rowAnnotation(Kingdom = taxonomy_core$Kingdom, 
                        col = list(Kingdom = color_mapping_rows),
                        gp = gpar(col = "black"),
                        simple_anno_size = unit(0.3, "cm"),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

### Creating heatmap
Heatmap(core_matrix, 
        show_column_names = FALSE, row_labels = core_species,
        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
        col = col_heatmap,
        name = "log10(1+x)",
        column_title = "Samples",
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 10),
        right_annotation = row_ha)

## Saving heatmap
png(file=paste0(path_plots, "/Micro_Heatmap.png"),
    width=30,height=40,units="cm",res=1200)
ht <- Heatmap(core_matrix, 
              show_column_names = FALSE, row_labels = core_species,
              row_names_gp = gpar(fontsize = 7, fontface = "italic"),
              col = col_heatmap,
              name = "log10",
              column_title = "Samples",
              column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 10),
              right_annotation = row_ha)
draw(ht, merge_legend = TRUE)
dev.off()
