##### !! Running SPIEC-EASI for co-occurrence analysis !! #####

## Load libraries
set.seed(1234)
library(phyloseq)
library(SpiecEasi)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(pals)


##### Input #####
pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"
cancer_type <- "TCGA-OV"

## Reading in manually reviewed core microorganisms
micro_fil <- readRDS(paste0(pathout, "/microorganisms_filtered_core.rds"))


##### Performing SPIEC-EASI with method glasso #####
####  (Tutorial from: https://biovcnet.github.io/_pages/NetworkScience_igraphviz.html )
spiec_glasso <- spiec.easi(micro_fil, method='glasso', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50, ncores = 5))

# Getting adjacency matrix and looking at numbers of edges
adj_mat <- getRefit(spiec_glasso)
table(as.numeric(adj_mat))

# Scaling the covariance matrix to become a correlation matrix and adding weights
corr_mat  <- cov2cor(as.matrix(getOptCov(spiec_glasso)))
weighted_adj <- corr_mat*getRefit(spiec_glasso)

## Creating the graph objects
grph_unweighted <- adj2igraph(adj_mat)
grph_weighted <- adj2igraph(weighted_adj)

# First look at the plots
plot(grph_unweighted, vertex.size=1, vertex.label=NA)
plot(grph_weighted, vertex.size=1, vertex.label=NA)

### Adding taxonomy data
taxonomy <- as.data.frame(tax_table(micro_fil))
V(grph_weighted)$Kingdom <- c(taxonomy$Kingdom)
V(grph_weighted)$Kingdom <- sub("Eukaryota", "Fungi", V(grph_weighted)$Kingdom)
V(grph_weighted)$Phylum <- c(taxonomy$Phylum)
V(grph_weighted)$Genus <- c(taxonomy$Genus)
V(grph_weighted)$Species <- c(taxonomy$Species)
V(grph_weighted)$name <- rownames(taxonomy)

### Removing nodes with no edges
rem_nodes_glasso <- names(which(igraph::degree(grph_weighted)<1))
cat("Removal of ", length(rem_nodes_glasso), "nodes, with no edges. Which are: \n")
print(rem_nodes_glasso)

## Saving csv of taxonomy of microbes with no interaction
rem_taxonomy_glasso <- taxonomy[rownames(taxonomy) %in% rem_nodes_glasso, ]
write.csv(rem_taxonomy_glasso, paste0(pathout, "/Removed_nodes_SPIECEASI_glasso.csv"))

## Removal
grph_weighted <- delete.vertices(grph_weighted, which(igraph::degree(grph_weighted)<1))

### Adjusting plot parameters
V(grph_weighted)$size <- (igraph::degree(grph_weighted) + 0.01) # Adjusting size based on number of edges
V(grph_weighted)$color <- ifelse(V(grph_weighted)$Kingdom == "Bacteria", "#377eb8",
                           ifelse(V(grph_weighted)$Kingdom == "Fungi", "#984ea3",
                                  ifelse(V(grph_weighted)$Kingdom == "Viruses", "#e41a1c", "black"))) # Setting color of nodes
E(grph_weighted)$color <- "#66c2a5" # Setting color of edges
E(grph_weighted)$color[E(grph_weighted)$weight<0] <- "#f46d43" # Setting color of negative edges
E(grph_weighted)$width <- abs(E(grph_weighted)$weight)*15 # Adjusting edge width based on weight

## Plotting the refined plot
png(filename = paste0(pathout, "/Plots/Glasso_Network.png"), width = 20, height = 20, units = "cm", res = 1200)

plot(grph_weighted,
     main= paste("Network of", cancer_type, sep = " "),
     layout=layout_with_dh(grph_weighted),
     vertex.label = NA
     #vertex.label = ifelse(V(grph_weighted)$size < 15, NA, V(grph_weighted)$Genus),
     #vertex.label.color = "black",
     #vertex.label.cex=0.8
     )
legend("topright", legend = c("Bacteria", "Fungi", "Viruses"), 
       col = c("#377eb8", "#984ea3", "#e41a1c"), pch = 16, cex = 0.8)

dev.off()

### Plotting the weighted matrix to get clusters
# Getting weighted adjacency matrix from filtered graph object
weighted_adjacency_matrix <- as.matrix(get.adjacency(grph_weighted, attr = "weight", sparse = FALSE))

# Getting taxonomy for kept species
kept_taxonomy_glasso <- taxonomy[! rownames(taxonomy) %in% rem_nodes_glasso, ]
kept_taxonomy_glasso$Kingdom <- sub("Eukaryota", "Fungi", kept_taxonomy_glasso$Kingdom)

## Row annotation
# Getting colors for each phyla 
color_mapping_rows <- setNames(c("#377eb8", "#984ea3", "#e41a1c"), c("Bacteria", "Fungi", "Viruses"))

# Creating row annotation
row_ha <- rowAnnotation(Kingdom = kept_taxonomy_glasso$Kingdom, 
                        col = list(Kingdom = color_mapping_rows),
                        gp = gpar(col = "black"),
                        simple_anno_size = unit(0.3, "cm"),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

## Column annotation
# Getting colors for each phyla 
unique_phyla <- unique(kept_taxonomy_glasso$Phylum)
num_phyla <- length(unique_phyla)

colors_columns <- trubetskoy(n = num_phyla)
color_mapping_columns <- setNames(colors_columns, unique_phyla)

# Creating column annotation
col_ha <- columnAnnotation(Phylum = kept_taxonomy_glasso$Phylum, 
                        col = list(Phylum = color_mapping_columns),
                        gp = gpar(col = "black"),
                        simple_anno_size = unit(0.3, "cm"),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

# Creating species names
species_labels = ifelse(kept_taxonomy_glasso$Kingdom == "Bacteria", paste(kept_taxonomy_glasso$Genus, kept_taxonomy_glasso$Species),
                      ifelse(kept_taxonomy_glasso$Kingdom == "Fungi", paste(kept_taxonomy_glasso$Genus, kept_taxonomy_glasso$Species),
                             ifelse(kept_taxonomy_glasso$Kingdom == "Viruses", paste(kept_taxonomy_glasso$Genus, kept_taxonomy_glasso$Species), "unknown kingdom")))

### Putting heatmap together
Heatmap(weighted_adjacency_matrix,
        name = "Weighted \nadjacency",
        col = colorRamp2(c(-0.3, 0, 0.3), c("turquoise", "black", "yellow")),
        right_annotation = row_ha,
        row_labels = species_labels,
        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
        bottom_annotation = col_ha,
        column_labels = species_labels,
        column_names_gp = gpar(fontsize = 7, fontface = "italic"))

## Saving heatmap
png(file=paste0(pathout, "/Plots/Glasso_Heatmap.png"),
    width=30,height=30,units="cm",res=1200)
ht <- Heatmap(weighted_adjacency_matrix,
              name = "Weighted \nadjacency",
              col = colorRamp2(c(-0.3, 0, 0.3), c("turquoise", "black", "yellow")),
              right_annotation = row_ha,
              row_labels = species_labels,
              row_names_gp = gpar(fontsize = 5, fontface = "italic"),
              bottom_annotation = col_ha,
              column_labels = species_labels,
              column_names_gp = gpar(fontsize = 5, fontface = "italic"),
              width = unit(10, "cm"),
              height = unit(10, "cm"))
draw(ht, merge_legend = TRUE)
dev.off()


#### trying to cluster (k-means)
library(cluster)
library(factoextra)

#create plot of number of clusters vs total within sum of squares
fviz_nbclust(weighted_adjacency_matrix, kmeans, method = "wss", print.summary = T)

# Needs changing according to elbow plot
number_of_clusters <- 8

# Plotting heatmap with clusters
group = kmeans(t(weighted_adjacency_matrix), centers = number_of_clusters, iter.max = 1000)$cluster
Heatmap(weighted_adjacency_matrix,
        cluster_columns = cluster_within_group(weighted_adjacency_matrix, group),
        cluster_rows = cluster_within_group(weighted_adjacency_matrix, group),
        name = "Weighted \nadjacency",
        col = colorRamp2(c(-0.3, 0, 0.3), c("red", "black", "cyan")),
        show_row_names = F,
        show_column_names = F,
        left_annotation = row_ha,
        top_annotation = col_ha)

png(file=paste0(pathout, "/Plots/Glasso_Heatmap_Cluster.png"),
    width=30,height=30,units="cm",res=1200)
ht_clust <- Heatmap(weighted_adjacency_matrix,
                    cluster_columns = cluster_within_group(weighted_adjacency_matrix, group),
                    cluster_rows = cluster_within_group(weighted_adjacency_matrix, group),
                    name = "Weighted \nadjacency",
                    col = colorRamp2(c(-0.3, 0, 0.3), c("red", "black", "cyan")),
                    show_row_names = F,
                    show_column_names = F,
                    left_annotation = row_ha,
                    top_annotation = col_ha)
draw(ht_clust, merge_legend = TRUE)
dev.off()

# Get groups as data frame
microbe_groups <- as.data.frame(t(data.frame(as.list(group))))
colnames(microbe_groups) <- "Group"
row.names(microbe_groups) <- sub("X", "", row.names(microbe_groups))

# Adding taxonomy information
microbe_groups <- merge(taxonomy, microbe_groups, by = "row.names")
colnames(microbe_groups)[1] <- "Tax_ID"

# Sorting based on group
microbe_groups <- microbe_groups[order(microbe_groups$Group), ]

# Save data frames with groups
write.csv(microbe_groups, paste0(pathout, "/SPIEC_groups.csv"), row.names = F)
