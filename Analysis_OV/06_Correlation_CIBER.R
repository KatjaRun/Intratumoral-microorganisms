##### !! Correlation between core microbiome and CIBERSORT !! #####

## Input
pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"

## Libraries
set.seed(1234)
library(CIBERSORT)
library(biomaRt)
library(phyloseq)
library(microbiome)
library(Hmisc)
library(tibble)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(pals)
library(RColorBrewer)
library(igraph)


##### [1] Performing CIBERSORT 
## Getting the LM22 signature matrix
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")

## Getting the raw counts
raw_counts <- read.csv(paste0(pathout, "/counts_raw_unstrand.csv"))

## Getting HUGO gene name from ensembl gene id
# Connecting to the Ensembl BioMart database
ensembl <- useMart("ensembl")

# Specify the dataset: human genes
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Getting Ensembl gene IDs and removing version information
raw_counts$X <- c(gsub("\\..*$", "", raw_counts$X))
ensembl_ids <- c(raw_counts$X)

# Specify the attributes you want to retrieve (Ensembl gene ID and HUGO gene name)
attributes <- c("ensembl_gene_id", "hgnc_symbol")

# Get the conversion
conversion <- getBM(attributes = attributes, filters = "ensembl_gene_id", 
                    values = ensembl_ids, mart = dataset)

# Merging the original data frame with the conversion data frame
raw_counts_HUGO <- merge(conversion, raw_counts, by.x = "ensembl_gene_id", by.y = "X", all.x = T)

# Getting it into the right format for CIBERSORT
names(raw_counts_HUGO)[names(raw_counts_HUGO) == "hgnc_symbol"] <- "Name"
raw_counts_HUGO <- raw_counts_HUGO[, -which(names(raw_counts_HUGO) == "ensembl_gene_id")]

# Filtering out LM22 genes to avoid problems with duplicate row names
data(LM22)
LM22_genes <- c(row.names(LM22))

raw_counts_LM22 <- raw_counts_HUGO[raw_counts_HUGO$Name %in% LM22_genes, ]

# Summing double rows
raw_counts_LM22_summed <- raw_counts_LM22 %>%
  group_by(Name) %>%
  summarise_all(sum)

# Saving the mixture file
mixture_file <- paste0(pathout, "/mixture_CIBERSORT.txt")
write.table(raw_counts_LM22_summed, mixture_file, sep = "\t", quote = FALSE, row.names = FALSE)


### Performing CIBERSORT
ciber <- cibersort(sig_matrix, mixture_file)

## Making the sample names match the ones from the phyloseq
row.names(ciber) <- substr(row.names(ciber), 1, 16)
row.names(ciber) <- gsub("\\.", "-", row.names(ciber))

## Saving the results of CIBERSORT
write.csv(ciber, paste0(pathout, "/Results_CIBERSORT.csv"))



#### [2] Performing the correlation
### Preparing the data frame for correlation
## Loading phyloseq of core microbiome
core <- readRDS(paste0(pathout, "/microorganisms_filtered_core.rds"))

## Transforming abundance data to clr
core_clr <- microbiome::transform(core, "clr")

## Getting OTU table with samples as rows
core_OTU <- as.data.frame(t(otu_table(core_clr)))

## Merging output from CIBERSORT and OTU table
data <- merge(core_OTU, ciber, by.x = "row.names", by.y = "X")

# Assigning samples as rownames
rownames(data) <- data$Row.names
data <- data[ , -1]

### Performing Pearson correlation with P-values
pear_corr <- rcorr(as.matrix(data))

### Getting parameters and adjusting the p-value
adj_matrix <- as.data.frame(pear_corr$r)

# Saving unfiltered correlation matrix
write.csv(adj_matrix, paste0(pathout, "/Correlation_Matrix.csv"))

# Adjusting the P value matrix
p_matrix <- as.data.frame(pear_corr$P)
p_adj_matrix <- matrix(p.adjust(as.vector(as.matrix(p_matrix)), method='fdr'),ncol=ncol(p_matrix))

# Adding back row and column names
colnames(p_adj_matrix) <- colnames(p_matrix)
rownames(p_adj_matrix) <- rownames(p_matrix)

##### [2.1] Drawing heatmap
# Filtering for significant P-values
cell_types <- colnames(ciber)
cell_types <- cell_types[-c(1, (length(cell_types) - 2):length(cell_types))]
microbes <- colnames(core_OTU)

corr_P_fil <- as.data.frame(p_adj_matrix[rownames(p_adj_matrix) %in% microbes, 
                                         colnames(p_adj_matrix) %in% cell_types])

P_ht <- corr_P_fil %>% rownames_to_column(var = "Tax_ID") %>%
  pivot_longer(cols = cell_types,
               names_to = "Immune_cells",
               values_to = "P_value")

P_sig <- P_ht[P_ht$P_value < 0.05, ]
sig_species <- unique(P_sig$Tax_ID)

## Filtering the correlation table for significant species
corr_r_ht <- adj_matrix[row.names(adj_matrix) %in% sig_species, colnames(adj_matrix) %in% cell_types]
corr_p_ht <- p_adj_matrix[row.names(adj_matrix) %in% sig_species, colnames(adj_matrix) %in% cell_types]

## Adjusting names before plotting
# Adding species names
core_tax <- as.data.frame(tax_table(core))
core_tax$Names <-  ifelse(core_tax$Kingdom %in% c("Bacteria", "Eukaryota"), 
                          paste(core_tax$Genus, core_tax$Species, sep = " "),
                          ifelse(core_tax$Kingdom == "Viruses",
                                 ifelse(nchar(core_tax$Species) == 1, 
                                        paste(core_tax$Genus, core_tax$Species, sep = " "), 
                                        core_tax$Species),
                                 "Unknown_Kingdom"))

filtered_core_tax <- core_tax[row.names(core_tax) %in% sig_species, ]

# Function for renaming
rename_species <- function(old_name, new_name) {
  filtered_core_tax$Names[filtered_core_tax$Names == old_name] <- new_name
  return(filtered_core_tax)
}

filtered_core_tax <- rename_species("haeme3", "Gemykibivirus haeme3")
filtered_core_tax <- rename_species("Pseudomonas uncultured Pseudomonas sp.","Pseudomonas sp.")
filtered_core_tax <- rename_species("Ambispora uncultured Glomeromycotina", "Glomeromycotina sp.")
filtered_core_tax <- rename_species("Severe acute respiratory syndrome-related coronavirus", "SARS-Cov")
filtered_core_tax <- rename_species("Moniliella uncultured Basidiomycota", "Moniliealla sp.")
filtered_core_tax <- rename_species("Desulfosudis Deltaproteobacteria bacterium", "Desulfosudis sp.")
filtered_core_tax <- rename_species("Lipomyces uncultured Ascomycota", "Lipomyces sp.")
filtered_core_tax <- rename_species("Oedogoniomyces uncultured Chytridiomycota", "Oedogoniomyces sp.")

species_names <- rownames(filtered_core_tax)
names(species_names) <- filtered_core_tax$Names

# Bringing filtered core tax in the right order before changing rownames
core_tax_ordered <- filtered_core_tax[match(rownames(corr_r_ht), rownames(filtered_core_tax)), ]
rownames(corr_r_ht) <- names(species_names)[match(rownames(corr_r_ht), as.character(species_names))]

## Changing names for immune cells
immune_names <- colnames(ciber)

immune_names <- immune_names[!(immune_names %in% c("X", "P.value", "Correlation", "RMSE"))]
names(immune_names) <- c("Naive B cells", "Memory B cells", "Plasma cells", "CD8 T cells", "Naive CD4 T", 
                         "Resting memory CD4 T", "Activated memory CD4 T", "Follicular helper T",
                         "Regulatory T cells", "Gamma delta T cells", "Resting NK cells", "Activated NK", 
                         "Monocytes", "M0 Macrophages", "M1 Macrophages", "M2 Macrophages", "Resting dendritic cells",
                         "Activated DC", "Resting mast cells", "Activated mast cells", "Eosinophils", 
                         "Neutrophils")

colnames(corr_r_ht) <- names(immune_names)[match(colnames(corr_r_ht), as.character(immune_names))]

## Creating row annotations
# For kingdoms
color_mapping_kingdoms <- setNames(c("#377eb8", "#984ea3", "#e41a1c"), c("Bacteria", "Eukaryota", "Viruses"))

col_Kingdom <- rowAnnotation(Kingdom = core_tax_ordered$Kingdom, 
                             col = list(Kingdom = color_mapping_kingdoms),
                             gp = gpar(col = "black"),
                             simple_anno_size = unit(0.3, "cm"),
                             annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                             annotation_legend_param = list(labels_gp = gpar(fontsize = 17)))

# For phyla
unique_phyla <- unique(core_tax_ordered$Phylum)
num_phyla <- length(unique_phyla)

colors_phyla <- trubetskoy(n = num_phyla)
color_mapping_phyla <- setNames(colors_phyla, unique_phyla)

col_Phyla <- rowAnnotation(Phylum = core_tax_ordered$Phylum, 
                           col = list(Phylum = color_mapping_phyla),
                           gp = gpar(col = "black"),
                           simple_anno_size = unit(0.3, "cm"),
                           annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                           annotation_legend_param = list(labels_gp = gpar(fontsize = 17)))


## Drawing and saving the heatmap
ht <- Heatmap(as.matrix(corr_r_ht),
              name = "Pearson \ncorrelation",
              col = colorRamp2(c(-0.5, -0.25, 0, 0.25, 0.5), c("#67001f", "#b2182b", "white", "#2166ac", "#053061")),
              row_names_gp = gpar(fontsize = 15, fontface = "italic"),
              column_names_gp = gpar(fontsize = 17),
              right_annotation = c(col_Kingdom, col_Phyla),
              column_km = 2, 
              row_km = 2,
              heatmap_legend_param = list(labels_gp = gpar(fontsize = 17)),
              cell_fun = function(j, i, x, y, w, h, fill) {
                asterisk_position_y <- y - h/2.5
                
                ifelse(corr_p_ht[i, j] < 0.001, grid.text("***", x, asterisk_position_y, gp = gpar(fontsize = 15)),
                       ifelse(corr_p_ht[i, j] < 0.01, grid.text("**", x, asterisk_position_y, gp = gpar(fontsize = 15)),
                              ifelse(corr_p_ht[i, j] < 0.05, grid.text("*", x, asterisk_position_y, gp = gpar(fontsize = 15)), 
                                     grid.text(""))))
              }
)

png(paste0(pathout, "/Plots/Immune_correlation_heatmap_NEW.png"), width=30,height=40,units="cm",res=1000)
ht_draw <- draw(ht, heatmap_legend_side = "left", merge_legend = T, padding = unit(c(10, 2, 2, 40), "mm"))
dev.off()

##### [2.2] Plotting the network
# Filtering adjacency matrix based on significant values
p_sig_matrix <- ifelse(p_adj_matrix < 0.05, TRUE, FALSE)

adj_network <- adj_matrix

adj_network[!p_sig_matrix] <- 0

## Replacing immune cell names in columns and rows
colnames(adj_network) <- ifelse(!(colnames(adj_network) %in% immune_names), colnames(adj_network),
                               names(immune_names)[match(colnames(adj_network), 
                                                         as.character(immune_names))])

row.names(adj_network) <- ifelse(!(row.names(adj_network) %in% immune_names), row.names(adj_network),
                                names(immune_names)[match(row.names(adj_network), 
                                                          as.character(immune_names))])

## Filtering out species with  no significant interaction with immune cells
adj_network <- adj_network[row.names(adj_network) %in% c(sig_species, names(immune_names)),
                         colnames(adj_network) %in% c(sig_species, names(immune_names))]

## Filtering out microbe-microbe correlations
adj_network[colnames(adj_network) %in% sig_species & row.names(adj_network) %in% sig_species, 
           colnames(adj_network) %in% sig_species & row.names(adj_network) %in% sig_species] <- 0

## Filtering out immune-immune correlations
adj_network[colnames(adj_network) %in% names(immune_names) & row.names(adj_network) %in% names(immune_names), 
           colnames(adj_network) %in% names(immune_names) & row.names(adj_network) %in% names(immune_names)] <- 0

## Saving filtered correlation matrix
write.csv(adj_network, paste0(pathout, "/Correlation_Matrix_fil.csv"))

## Creating graph object
network <- igraph::graph_from_adjacency_matrix(as.matrix(adj_network), mode="undirected", diag = FALSE, 
                                               weighted = TRUE)

### Adding taxonomy data
cols_network <- colnames(adj_network)
cols_info <- data.frame("Info" = cols_network,
                        "Origin" = NA)
cols_info$Origin <- ifelse(cols_info$Info %in% names(immune_names), "Immune cell", "Microbe")
cols_info_sub <- cols_info[cols_info$Origin == "Microbe", ]

cols_info_sub <- merge(cols_info_sub, filtered_core_tax, by.x = "Info", by.y = "row.names", sort = F)

cols_info_microbes <- data.frame("Info" = cols_info_sub$Names,
                                 "Origin" = cols_info_sub$Kingdom)

cols_info_sub2 <- cols_info[cols_info$Origin == "Immune cell", ]

cols_info <- rbind(cols_info_microbes, cols_info_sub2)

V(network)$Kingdom <- c(cols_info$Origin)
V(network)$Kingdom <- sub("Eukaryota", "Fungi", V(network)$Kingdom)
V(network)$Species <- c(cols_info$Info)

### Adjusting plot parameters
V(network)$size <- (igraph::degree(network)) # Adjusting size based on number of edges
V(network)$color <- ifelse(V(network)$Kingdom == "Bacteria", "#377eb8",
                           ifelse(V(network)$Kingdom == "Fungi", "#984ea3",
                                  ifelse(V(network)$Kingdom == "Viruses", "#e41a1c", 
                                         ifelse(V(network)$Kingdom == "Immune cell", "seagreen",
                                                "black")))) # Setting color of nodes
E(network)$color <- "skyblue" # Setting color of edges
E(network)$color[E(network)$weight<0] <- "tomato" # Setting color of negative edges
E(network)$width <- abs(E(network)$weight)*10 # Adjusting edge width based on weight

## Plotting the refined plot
layout <- layout_with_dh(network)
x_coordinates <- layout[, 1]
y_coordinates <- layout[, 2]

new_x_coordinates <- x_coordinates #*30
new_y_coordinates <- y_coordinates #*30

V(network)$x <- new_x_coordinates
V(network)$y <- new_y_coordinates

png(filename = paste0(pathout, "/Plots/Correlation_Network_NEW.png"), width = 20, height = 20, units = "cm", res = 1200)

plot(network,
     main= "Network of significant immune/microbe pearson correlations",
     vertex.label = ifelse((V(network)$Kingdom == "Immune cell" & V(network)$size > 10), 
                           V(network)$Species, NA),
     vertex.label.color = "black",
     vertex.label.cex=1,
     vertex.label.font = 2) # 2 stands for bold

legend("topright", 
       legend = c("Bacteria", "Fungi", "Viruses", "Immune cells", "Positive Correlation", "Negative Correlation"), 
       col = c("#377eb8", "#984ea3", "#e41a1c", "seagreen", "skyblue", "tomato"), 
       pch = c(16, 16, 16, 16, 15, 15),  # Use 16 for points, 95 for square symbol
       cex = 0.8)

dev.off()


## Drawing plot for subset of immune cells
# Filtering network to only nodes connected to three immune cells with highest number of edges
# Which are: M2 Macrophages, DCs, Memory B cells
relevant_nodes <- c("M2 Macrophages", "Activated NK", "Memory B cells", "Activated DC")
relevant_nodes <- union(relevant_nodes, neighbors(network, "M2 Macrophages")$name)
relevant_nodes <- union(relevant_nodes, neighbors(network, "Activated NK")$name)
relevant_nodes <- union(relevant_nodes, neighbors(network, "Memory B cells")$name)
relevant_nodes <- union(relevant_nodes, neighbors(network, "Activated DC")$name)
relevant_network <- induced_subgraph(network, relevant_nodes)

# Filtering out species which are only connected to one of those immune cells
rem_nodes_subset <- names(which(igraph::degree(relevant_network) == 1))
cat("Removal of ", length(rem_nodes_subset), "nodes, with no edges. Which are: \n")
print(rem_nodes_subset)

# Removal
relevant_network <- delete.vertices(relevant_network, which(igraph::degree(relevant_network) == 1))

# Plot the filtered network
png(filename = paste0(pathout, "/Plots/Correlation_Network_Subs.png"), width = 20, height = 20, units = "cm", 
    res = 1200)

plot(relevant_network,
     main = "Subset of significant immune/microbe pearson correlations",
     vertex.label = ifelse(V(relevant_network)$Kingdom == "Immune cell", 
                           V(relevant_network)$Species, 
                           ifelse((V(relevant_network)$Kingdom != "Immune cell" & V(relevant_network)$size > 5), 
                                  V(relevant_network)$Species, NA)),
     vertex.label.color = "black",
     vertex.label.cex = ifelse(V(relevant_network)$Kingdom == "Immune cell", 0.8, 0.8),
     vertex.label.font = ifelse(V(relevant_network)$Kingdom == "Immune cell", 2, 3)) # bold for immune cells and italic for species

legend("topleft", 
       legend = c("Bacteria", "Fungi", "Viruses", "Immune cells", "Positive Correlation", "Negative Correlation"), 
       col = c("#377eb8", "#984ea3", "#e41a1c", "seagreen", "skyblue", "tomato"), 
       pch = c(16, 16, 16, 16, 15, 15),  # Use 16 for points, 95 for square symbol
       cex = 0.8)

dev.off()









