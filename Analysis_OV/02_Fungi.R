#### !! Pipeline for fungi !! ####

#### Libraries
set.seed(43)
library(phyloseq)
library(speedyseq)
library(decontam)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(pals)
library(ComplexHeatmap)

#### Input 
pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"
path_fungi <- paste0(pathout, "/fungi.rds")
path_plots <- paste0(pathout, "/Plots")


##### [1] Load Phyloseq object #####
fungi <- readRDS(path_fungi)

## Creating a table with easy accessible species names
taxdata <- as.data.frame(tax_table(fungi))
species_names <- data.frame(tax_id = row.names(taxdata),
                            species_full = paste(taxdata$Genus, taxdata$Species, sep = " "),
                            species_short = paste(substr(taxdata$Genus, 1, 1), ". ", taxdata$Species, sep = ""))


##### [2] Decontamination #####
### [2.1] Filtering out unclassified genera and singletons
## Filtering out unclassified genera
fungi_decontam <- subset_taxa(fungi, Genus != "")

## Filtering out Singletons (taxa with only one read)
if (sum(taxa_sums(fungi_decontam) == 1) == 0) {
  print("No taxa with count equal to 1. Skipping prune_taxa.")
} else {
  cat("Pruning taxa with count equal to 1: ", sum(taxa_sums(fungi_decontam) == 1))
  fungi_decontam <- prune_taxa(taxa_sums(fungi_decontam) >1, fungi_decontam)
}


# Calculating how many taxa were lost
cat("Taxa lost after singleton and unclassified genera removal:", 
    1 - ntaxa(fungi_decontam)/ntaxa(fungi))


### [2.2] Filtering by Poore et al Plate-Center
## Function for extracting last n characters from R string
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

## Getting information about aliquots from sample_data
sampledata <- as.data.frame(sample_data(fungi_decontam))

# Defining Plate-Center combinations with at least 10 samples
tmp <- as.character(sampledata$bcr_aliquot_barcode)
sampledata$PlateCenter <- factor(substrRight(tmp, 7))

booleanPlateCenter <- as.logical(table(sampledata$PlateCenter) >10)
sufficientPlateCenter <- names(table(sampledata$PlateCenter))[booleanPlateCenter]

# Filtering out the sufficient combinations
sampledata$PlateCenterFlag <- (sampledata$PlateCenter %in% sufficientPlateCenter)
sampledata_PlateCenterSubset <- subset(sampledata, PlateCenterFlag)

# Creating a filtered phyloseq object
fungi_PlateCenterSubset <- subset_samples(fungi_decontam, sample_names(fungi_decontam) %in% 
                                               sampledata_PlateCenterSubset$sample_id)

cat("Insufficient samples:", nsamples(fungi_decontam)-nsamples(fungi_PlateCenterSubset))

# Applying decontam to filtered Phyloseq
decontam_fungi <- isContaminant(seqtab = fungi_PlateCenterSubset,
                                   conc = sampledata_PlateCenterSubset$concentration,
                                   method = "frequency",
                                   batch = sampledata_PlateCenterSubset$PlateCenter,
                                   threshold = 0.1)

# Printing out contaminated taxa
table(decontam_fungi$contaminant)

# Pruning taxa based on the identified contaminants
fungi_decontam <- prune_taxa(!decontam_fungi$contaminant, fungi_decontam)

# Calculating how much data was lost since the beginning
cat("Lost taxa after decontam by Plate-Center: ", 1 - ntaxa(fungi_decontam)/ntaxa(fungi))

## Saving table of removed taxa
taxa_contam <- decontam_fungi[decontam_fungi$contaminant == "TRUE", ]
taxa_contam$tax_id <- row.names(taxa_contam)
taxa_contam <- merge(species_names, taxa_contam, by = "tax_id")

write.csv(taxa_contam, file = paste(pathout, "/removed_contam_fungi.csv", sep = ""), row.names = F)

##### Saving the decontaminated output #####
saveRDS(fungi_decontam, file = paste(pathout, "/fungi_decontam.rds", sep = ""))

cat("Number of taxa at beginning: ", ntaxa(fungi))
cat("Number of taxa after decontamination: ", ntaxa(fungi_decontam))
cat("Number of taxa lost: ", ntaxa(fungi) - ntaxa(fungi_decontam))


##### [3] Alpha diversity #####

### Creating table with diversity parameters
## Should be used with the count data
# Calculating diversity, richness and rarity of the samples
diversity <- microbiome::diversity(fungi_decontam)
richness <- microbiome::richness(fungi_decontam)
rarity <- microbiome::rarity(fungi_decontam)
dominance <- microbiome::dominance(fungi_decontam)

# Combining all data frames and adding a constant x-axis for plotting
fungi_alpha <- cbind(diversity, richness, rarity, dominance)
fungi_alpha$constant_x <- 0

## Creating ViolinPLot for observed species richness
observed <- ggplot(fungi_alpha, aes(x = constant_x, y = observed)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Observed species richness") +
  ylim(0, 10000) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for species richness (Chao1)
chao1 <- ggplot(fungi_alpha, aes(x = constant_x, y = chao1)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("chao1 species richness") +
  ylim(0, 10000) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for Shannon diversity
shannon <- ggplot(fungi_alpha, aes(x = constant_x, y = shannon)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Shannon diversity") +
  ylim(0, 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for inverse simpson
inv_simpson <- ggplot(fungi_alpha, aes(x = constant_x, y = inverse_simpson)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Inverse Simpson Index") +
  ylim(0, 100) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for dominance core abundance
dominance_core <- ggplot(fungi_alpha, aes(x = constant_x, y = core_abundance * 100)) +
  geom_violin(trim = FALSE, fill = "skyblue") +
  geom_boxplot(width = 0.05, fill = "white") +
  # geom_jitter(color = "midnightblue", alpha = 0.5) +
  theme_light() +
  xlab("Samples") +
  ylab("Dominance of core abundance [%]") +
  ylim(0, 100) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Creating ViolinPlot for low abundance
low_abundance <- ggplot(fungi_alpha, aes(x = constant_x, y = low_abundance * 100)) +
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

ggsave(paste(path_plots, "/Fungi_alpha.png", sep = ""),
       width = 22, height = 25.5, units = "cm")



##### [4] Beta diversity #####

### Ordination using NMDS and bray distance
ord_NMDs <- phyloseq::ordinate(fungi_decontam, "NMDS", "bray")

options(repr.plot.width = 3, repr.plot.height = 3)
plot_ordination(fungi_decontam, ord_NMDs, type = "samples", color = "Subtype_mRNA", shape = "sample_type",
                title = "Fungi: NMDs of samples") +
  theme_light() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1)

ggsave(paste0(path_plots, "/Fungi_NMDS.png"), height = 10, width = 15, unit = "cm")

### Ordination using PCoA and bray distance
ord_PCOA <- ordinate(fungi_decontam, "PCoA")

plot_ordination(fungi_decontam, ord_PCOA, type = "samples", color = "Subtype_mRNA", shape = "sample_type",
                title = "Fungi: PCoA of samples") +
  theme_light() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1)

ggsave(paste0(path_plots, "/Fungi_PCoA.png"), height = 10, width = 15, unit = "cm")



##### [5] Abundance plots #####
## Calculating relative abundance
fungi_rel = microbiome::transform(fungi_decontam, "compositional")

## Plot for top 20 phyla
# Performing tax_glom at phylum level and saving 
fungi_phy = speedyseq::tax_glom(fungi_rel, taxrank = "Phylum")
saveRDS(fungi_phy, paste0(pathout, "/fungi_phyla.rds"))

# Sorting to get most abundant phyla and pruning the phyloseq
top20_phy = names(sort(taxa_sums(fungi_phy), decreasing = TRUE))[1:20]
fungi_top20_phy = prune_taxa(top20_phy, fungi_phy)

# Checking relative abundance of phyla
cat("Mean abundance of phyla: ", mean(colSums(otu_table(fungi_top20_phy))))
cat("Ranging from: ", min(colSums(otu_table(fungi_top20_phy))), "to", 
    max(colSums(otu_table(fungi_top20_phy))))

# Preparing data frame for plotting
top_phyla_OTU <- as.data.frame(otu_table(fungi_top20_phy))
top_phyla_tax <- as.data.frame(tax_table(fungi_top20_phy))

phylum_rownames <- c(top_phyla_tax$Phylum)
top_phyla_OTU$phyla <- phylum_rownames

tidy_phyla <- pivot_longer(top_phyla_OTU, cols = -phyla, names_to = "Sample", values_to = "Abundance")

# Preparing colors for plot
# Getting colors from pals package
colors <- pals::trubetskoy(n = 20)

# Extracting only the hexadecimal color codes
colors <- rev(rgb(t(col2rgb(colors)), maxColorValue = 255))

# Creating and saving plot for top 20 phyla
ggplot(tidy_phyla, aes(x = reorder(Sample, Abundance), y = Abundance, fill = reorder(phyla, Abundance))) +
  geom_bar(stat="identity", position="stack") +
  xlab("Sample") +
  ylab("Relative abundance") +
  ggtitle("Fungi: Top 20 Phyla") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm")) + 
  theme(legend.text = element_text(size = 7)) +
  scale_fill_manual(values = c(colors), guide = guide_legend(reverse = TRUE))

ggsave(paste0(path_plots, "/Fungi_top_phyla.png"))  


## Plot for top 20 genera
# Performing tax_glom at genus level and saving
fungi_gen = speedyseq::tax_glom(fungi_rel, taxrank = "Genus")
saveRDS(fungi_gen, paste0(pathout, "/fungi_genera.rds"))

# Sorting to get most abundant phyla and pruning the phyloseq
top20_gen = names(sort(taxa_sums(fungi_gen), decreasing = TRUE))[1:20]
fungi_top20_gen = prune_taxa(top20_gen, fungi_gen)

# Checking the relative abundance per sample
cat("Mean abundance of genera: ", mean(colSums(otu_table(fungi_top20_gen))))
cat("Ranging from: ", min(colSums(otu_table(fungi_top20_gen))), "to", 
    max(colSums(otu_table(fungi_top20_gen))))

# Preparing data frame for plotting
top_genera_OTU <- as.data.frame(otu_table(fungi_top20_gen))
top_genera_tax <- as.data.frame(tax_table(fungi_top20_gen))

genera_rownames <- c(top_genera_tax$Genus)
top_genera_OTU$genera <- genera_rownames

tidy_genera <- pivot_longer(top_genera_OTU, cols = -genera, names_to = "Sample", values_to = "Abundance")

# Creating and saving plot for top 20 genera
ggplot(tidy_genera, aes(x = reorder(Sample, Abundance), y = Abundance, fill = reorder(genera, Abundance))) +
  geom_bar(stat="identity", position="stack") +
  xlab("Sample") +
  ylab("Relative abundance") +
  ggtitle("Fungi: Top 20 Genera") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm")) + 
  theme(legend.text = element_text(size = 7)) +
  scale_fill_manual(values = c(colors), guide = guide_legend(reverse = TRUE))

ggsave(paste0(path_plots, "/Fungi_top_genera.png"))


##### [6] Defining core #####

### Defining core fungil species (>0.0001 relative abundance and present in 80% of all samples)
## Tax glom at species level
fungi_spe = speedyseq::tax_glom(fungi_decontam, taxrank = "Species")

saveRDS(fungi_spe, paste0(pathout, "/fungi_species.rds"))

## Transforming data into relative abundance and setting parameters for filtering
fungi_spe_rel = microbiome::transform(fungi_spe, "compositional")
prev = 0.8
abund = (1/ntaxa(fungi_spe)) #* 2)
  
# Getting mean species abundance
total_abundance <- rowSums(otu_table(fungi_spe_rel), na.rm = TRUE)
cat("Mean total abundance: ", mean(total_abundance, na.rm = TRUE))
cat("Median total abundance: ", median(total_abundance, na.rm = TRUE))

## Defining the Core species
fungi_core = core(fungi_spe_rel, detection = abund, prevalence = prev, include.lowest = T)

## Getting species names and ids for heatmap
taxonomy_core <- as.data.frame(tax_table(fungi_core))
core_id <- row.names(taxonomy_core)
core_species = c(paste(taxonomy_core$Genus, taxonomy_core$Species))

## Getting sample data for molecular subtypes
sample_core <- as.data.frame(sample_data(fungi_core))

## Preparing matrix for heatmap
fungi_core_prune = prune_taxa(core_id, fungi_spe)

# Saving untransformed data of core species
saveRDS(fungi_core_prune, paste0(pathout, "/fungi_core.rds"))

## Creating log-transformed matrix for heatmap
fungi_core_log = microbiome::transform(fungi_core_prune, "log10")
core_matrix <- data.matrix(as.data.frame(otu_table(fungi_core_log)))

## Preparing colors for heatmap
col_heatmap <- rev(cividis(n = 255))

## Row annotation
# Getting colors for each phyla 
unique_phyla <- unique(taxonomy_core$Phylum)
num_phyla <- length(unique_phyla)

colors_rows <- cubicl(n = num_phyla)
color_mapping_rows <- setNames(colors_rows, unique_phyla)

# Creating row annotation
row_ha <- rowAnnotation(Phylum = taxonomy_core$Phylum, 
                        col = list(Phylum = color_mapping_rows),
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
png(file=paste0(path_plots, "/Fungi_Heatmap.png"),
    width=32,height=20,units="cm",res=1200)
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

## Creating a csv with all core species that can be manually evaluated
core_tax <- as.data.frame(tax_table(fungi_core_prune))

write.csv(core_tax, paste0(pathout, "/Fungi_core.csv"), row.names = T)


### Table of core species is getting manually reviewed by:
## 1) Checking presence in MiMeDB
## 2) Checking available literature