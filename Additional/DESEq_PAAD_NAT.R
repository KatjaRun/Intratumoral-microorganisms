##### DESEq PAAD #####
set.seed(1234)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(ggplot2)
library(ggpubr)

pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_Data/"
pathbiom <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/PNAT_Results/Final_nt.biom"

##### Loading BIOM file and adding sample data
data <- import_biom(pathbiom)

## Changing names of tax_table and removing prefixes
colnames(tax_table(data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
tax_table(data)[, colnames(tax_table(data))] <- gsub(tax_table(data)[, colnames(tax_table(data))], 
                                                     pattern = "[a-z]__", replacement = "")


## Only bacteria
bacteria = subset_taxa(data, Kingdom == "Bacteria")

## Only fungi
fungi = subset_taxa(data, grepl("mycota$|Microsporidia|Nephridiophagidae", Phylum))

## Only viruses
virus = subset_taxa(data, Kingdom == "Viruses")

## Putting together all microorganisms
microorganisms = merge_phyloseq(bacteria, fungi, virus)

# Removing taxa with no genus
micro_decontam <- subset_taxa(microorganisms, Genus != "")

# Creating sample data with tissue type
# Gsub _bracken_species from samples
sample_names(micro_decontam) <- gsub("_bracken_species", "", sample_names(micro_decontam))

# Getting sample names
samples <- sample_names(micro_decontam)

# Adding information about origin
sample_types <- ifelse(grepl("^.{13}11", samples), "Normal",
                       ifelse(grepl("^.{13}01", samples), "Tumor",
                              ifelse(grepl("^.{13}06", samples), "Metastasis", NA)))

# Creating sample_data and adding to phyloseq
sam <- sample_data(data.frame(row.names = samples, Type = sample_types, Pairs = substr(samples, 1, 12)))
micro_decontam <- merge_phyloseq(micro_decontam, sam)

## Filtering out Singletons (taxa with only one read)
if (sum(taxa_sums(micro_decontam) == 1) == 0) {
  print("No taxa with count equal to 1. Skipping prune_taxa.")
} else {
  cat("Pruning taxa with count equal to 1: ", sum(taxa_sums(micro_decontam) == 1))
  micro_decontam <- prune_taxa(taxa_sums(micro_decontam) >1, micro_decontam)
}

# Calculating how many taxa were lost
cat("Taxa lost after singleton and unclassified genera removal:", 
    1 - ntaxa(micro_decontam)/ntaxa(bacteria))


##### [1] DESEq for normal-tumor #####
### Filtering out the recurrent samples through the Sample Type Code
normal <- samples[grep("^.{13}11", samples)]
norm_patients <- substr(normal, 1, 12)

### Getting the primary samples for the recurrent samples
tumor_norm <- c()
for (i in norm_patients) {
  pattern <- paste0("^", i, "-01")
  tumor_norm <- c(tumor_norm, samples[grep(pattern, samples)])
}

### Filtering phyloseq for those samples
norm_micro <- subset_samples(micro_decontam, sample_names(micro_decontam) %in% c(tumor_norm, normal))

### Creating a DESeq2 object
norm_genus <- speedyseq::tax_glom(norm_micro, "Genus")

norm_micro_deseq <- phyloseq_to_deseq2(norm_genus, ~ Type)

### Performing DESeq
data_DESeq <- DESeq(norm_micro_deseq)

### Extracting the results
res <- results(data_DESeq)
summary(res)
res <- res[order(res$padj),]

### Getting significantly abundant species
alpha = 0.5
table(res$padj < alpha)

significant <- subset(res, padj < 0.5)

## Getting the species names from tax ID
sig_species <- row.names(significant)
taxonomy <- as.data.frame(tax_table(norm_genus))
sig_taxonomy <- taxonomy[rownames(taxonomy) %in% sig_species, ]

cat("Significant species with padj <0.5: ", paste(sig_taxonomy$Genus, sig_taxonomy$Species, sep = "_"))

### Generating a Volcano plot
par(mfrow=c(1,1))

png(filename = paste0(pathout, "/Volcano_Tum_Norm.png"), width = 15, height = 10, unit = "cm", res = 800)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
with(subset(res, padj<0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<0.5 & abs(log2FoldChange)>10), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

## CLR transformation of norm_genus
norm_clr <- microbiome::transform(norm_genus, "clr")

### Creating plot for Mycobacterium
myco <- data.frame(Samples = sample_names(norm_clr),
                   Pairs = substr(sample_names(norm_clr), 1, 12),
                   Category = sample_data(norm_clr)$Type,
                   Myco = t(otu_table(norm_clr)["1773",]))

## Reading in microorganisms of PAAD
PAAD <- readRDS("/data/projects/2020/OvarianCancerHH/Thesis_Katja/PAAD_Data_nt/microorganisms.rds")
PAAD_genus <- speedyseq::tax_glom(PAAD, "Genus")
PAAD_clr <- microbiome::transform(PAAD_genus, "clr")

PAAD_otu <- as.data.frame(otu_table(PAAD_clr))
PAAD_myco <- as.data.frame(t(PAAD_otu["1773", ]))
PAAD_myco$constant <- 0

# Plotting
comparison <- ggplot(myco, aes(Category, X1773, fill = Category)) + 
  geom_boxplot() +
  geom_point() + 
  geom_line(aes(group = Pairs)) +
  theme_bw() +
  ylab("Mycobacterium abundance (clr)") + xlab("") +
  scale_fill_manual(values = c("skyblue", "tomato")) +
  ylim(-0.5, 13) +
  theme(axis.text.x = element_text(size = 10, color = "black"),  
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())  +
  labs(fill = "Sample Type") +
  theme(plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  stat_compare_means(label = "p.signif")


all <- ggplot(PAAD_myco, aes(constant, `1773`)) + 
  geom_boxplot(fill = "tomato") +
  theme_bw() +
  ylab("") + xlab("") +
  ylim(-0.5, 13) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(breaks = 0, labels = "All Tumor") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"))


ggarrange(comparison, all, nrow = 1, common.legend = T, legend = "right", 
          widths = c(1, 0.5), align = "h")

ggsave(paste0(pathout, "/Mycobacterium_Norm_Tum.png"), width = 5.6, height = 3.2)


#### [1.2] Ordination plot for tumor-normal
### Ordination using PCoA and bray distance
ord_PCOA <- ordinate(norm_micro, "PCoA")

options(repr.plot.width = 3, repr.plot.height = 3)
plot_ordination(norm_micro, ord_PCOA, color = "Type", shape = "Pairs") +
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("Normal" = "skyblue", "Tumor" = "tomato"), name = "Sample Type") +
  scale_shape_manual(values = c(19, 17, 15, 18), name = "Patient") +
  geom_point(size = 4, alpha = 0.7)

ggsave(paste0(pathout, "/Beta_PCoA.png"), height = 10, width = 15, unit = "cm")



##### [2] DESEq for tumor-metastasis #####
### Filtering out the recurrent samples through the Sample Type Code
meta <- samples[grep("^.{13}06", samples)]
meta_patients <- substr(meta, 1, 12)

### Getting the primary samples for the recurrent samples
tumor_meta <- c()
for (i in meta_patients) {
  pattern <- paste0("^", i, "-01")
  tumor_meta <- c(tumor_meta, samples[grep(pattern, samples)])
}

### Filtering phyloseq for those samples
meta_micro <- subset_samples(micro_decontam, sample_names(micro_decontam) %in% c(tumor_meta, meta))


####Performing differential abundance analysis #####
### Creating a DESeq2 object
meta_genus <- speedyseq::tax_glom(meta_micro, "Genus")

meta_micro_deseq <- phyloseq_to_deseq2(meta_genus, ~ Type)

### Performing DESeq
data_DESeq <- DESeq(meta_micro_deseq)

### Extracting the results
res <- results(data_DESeq)
summary(res)
res <- res[order(res$padj),]

### Getting significantly abundant species
alpha = 0.5
table(res$padj < alpha)

significant <- subset(res, padj < 0.5)

## Getting the species names from tax ID
sig_species <- row.names(significant)
taxonomy <- as.data.frame(tax_table(meta_genus))
sig_taxonomy <- taxonomy[rownames(taxonomy) %in% sig_species, ]

cat("Significant species with padj <0.5: ", paste(sig_taxonomy$Genus, sig_taxonomy$Species, sep = "_"))

### Generating a Volcano plot
par(mfrow=c(1,1))

png(filename = paste0(pathout, "/Volcano_Tum_Meta.png"), width = 10, height = 10, unit = "cm", res = 800)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
with(subset(res, padj<0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<0.5 & abs(log2FoldChange)>10), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()




