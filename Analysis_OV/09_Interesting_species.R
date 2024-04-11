##### !! Filtering for most interesting species !! #####

pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"

## Libraries
library(phyloseq)
library(venneuler)
library(eulerr)


### [1] Loading species information
## Total core
core <- readRDS(paste0(pathout, "/microorganisms_filtered_core.rds"))
core_tax <- as.data.frame(tax_table(core))
core_tax$Names <-  ifelse(core_tax$Kingdom %in% c("Bacteria", "Eukaryota"), 
                          paste(core_tax$Genus, core_tax$Species, sep = " "),
                          ifelse(core_tax$Kingdom == "Viruses",
                                 ifelse(nchar(core_tax$Species) == 1, 
                                        paste(core_tax$Genus, core_tax$Species, sep = " "), 
                                        core_tax$Species),
                                 "Unknown_Kingdom"))
core_species <- core_tax$Names

## Correlation immune cells
immune <- read.csv(paste0(pathout, "/Correlation_Matrix_fil.csv"), row.names = 1)
immune_ids <- rownames(immune)[grepl("^\\d+$", rownames(immune))]
immune_tax <- core_tax[row.names(core_tax) %in% immune_ids, ]
immune_species <- immune_tax$Names

## Survival Analysis
surv <- read.table(file.path(pathout, "/Survival_Species.txt"), header = TRUE, sep = "\t")
surv_species_fil <- surv[surv$P < 0.05, ]
surv_species <- surv_species_fil$Species

## GSVA correlation
GSVA <- read.csv(paste0(pathout, "/GSVA_gene_sets/Regression_table_fil.csv"))
GSVA_species <- unique(GSVA$Species)
GSVA_clean <- na.omit(GSVA_species)
GSVA_clean <- as.vector(GSVA_clean)


##### [2] Drawing a venneuler diagram #####
## Creating a dataframe for venn diagram
species_from <- c("core_species", "immune_species", "surv_species", "GSVA_clean")
#species_from <- c("immune_species", "surv_species", "GSVA_species")
dfs <- list()

# Loop over each group
for (i in species_from) {
  spe <- get(i)  
  out <- data.frame(Species = spe, Group = i)  
  dfs[[length(dfs) + 1]] <- out 
}

# Combining all data frames 
venn_data <- bind_rows(dfs)

# Getting species combinations
combinations <- paste(venn_data$Species, venn_data$Group, sep = "+")

# Creating empty dataframe for T/F matrix
eulerr_df <- data.frame(matrix(nrow = length(core_species), ncol = length(species_from)))
rownames(eulerr_df) <- core_species
colnames(eulerr_df) <- species_from

# Getting T/F values
for (comb in combinations) {
  parts <- strsplit(comb, "\\+")
  species_name <- parts[[1]][1]
  group <- parts[[1]][2]
  eulerr_df[species_name, group] <- TRUE
}

# Replacing NA with FALSE
eulerr_df[is.na(eulerr_df)] <- FALSE

# Get better column names
colnames(eulerr_df) <- c("Core", "Cibersort", "Survival", "GSVA")
colnames(eulerr_df) <- c("Core", "Immune cells", "Survival", "Immune parameters")

# Performing euler fitting
fit <- euler(eulerr_df)
fit$stress
fit$diagError

# Plotting the eulerr diagram
png(filename = paste0(pathout, "/Plots/Eulerr_sigspecies.png"), width = 15, height = 7, units = "cm", res = 1200)
plot(fit, fill = c("snow", "gold", "seagreen", "darkorange1"), quantities = T, labels = F, legend = T)
dev.off()

### Getting overlapping bacteria
surv_immune <- intersect(surv_species, immune_species)
immune_GSVA <- intersect(immune_species, GSVA_species)

### Saving interesting species
SI_tax <- core_tax[core_tax$Names %in% surv_immune, ]
SI_tax$group <- "Surv_Immune"

IG_tax <- core_tax[core_tax$Names %in% immune_GSVA, ]
IG_tax$group <- "Immune_GSVA"

# Binding together
sig_tax <- rbind(SI_tax, IG_tax)

# Saving csv
write.csv(sig_tax, paste0(pathout, "/Important_species.csv"))


#### [2] Species for epitope prediction
app <- GSVA[GSVA$Gene_set == "GOBP_APP" & GSVA$Corr_r > 0, ]
cyt <- GSVA[GSVA$Gene_set == "CYT" & GSVA$Corr_r > 0, ]

app_species <- unique(app$Species)
cyt_species <- unique(cyt$Species)

