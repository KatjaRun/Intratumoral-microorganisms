#### !! Survival analysis !! ####

pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"

## Libraries
library(phyloseq)
library(compositions)
library(survival)
library(ggsurvfit)


### [1] Loading needed data
# 1. Survival data from UCSC XENA (from Script 01)
survival <- read.csv(paste0(pathout, "/Survival_data.csv"))

# 2. Clinical data
# Download process from GDC:
# Data category: Clinical, Data Format: bcr biotab
# ...clinical_patient_ov.txt -> Put download folder into pathout, extract file to be in that folder
clinical_folder <- list.files(path = pathout, pattern = "^gdc")

clinical_path <- file.path(pathout, clinical_folder, 
                           (list.files(path = file.path(pathout, clinical_folder), pattern = "^nationwidechildrens")))

clinical <- read.delim(file = clinical_path, sep = "\t")

# 3. CIBERSORT data
ciber <- read.csv(paste0(pathout, "/Results_CIBERSORT.csv")) 
ciber_perc <- subset(ciber, select = - c(X, P.value, Correlation, RMSE))
row.names(ciber_perc) <- ciber$X
ciber_perc <- ciber_perc * 100

# 4. Filtered core microbiome
core <- readRDS(file.path(pathout, "microorganisms_filtered_core.rds"))

# Performing clr normalization
core_clr <- microbiome::transform(core, "clr")


### [2] Merging data together
## Merging OTU and tax table
core_otu <- as.data.frame(otu_table(core_clr))
core_tax <- as.data.frame(tax_table(core_clr))
core_tax$Names <-  ifelse(core_tax$Kingdom %in% c("Bacteria", "Eukaryota"), 
                          paste(core_tax$Genus, core_tax$Species, sep = " "),
                          ifelse(core_tax$Kingdom == "Viruses",
                                 ifelse(nchar(core_tax$Species) == 1, 
                                        paste(core_tax$Genus, core_tax$Species, sep = " "), 
                                        core_tax$Species),
                                 "Unknown_Kingdom"))

# Giving meaningful rownames
rownames(core_otu) <- core_tax$Names

# Exchanging rows and columns for merging with survival data
core_otu <- as.data.frame(t(core_otu))

## Merging with survival data
# Stripping last character from sample names to get same length
core_otu$sample <- sub(".$", "", rownames(core_otu))
core_otu$Original_ID <- row.names(core_otu)

# Merging dataframes by samples 
core_survival <- merge(survival, core_otu, by = "sample")

# Subsetting the clinical dataframe
clinical_sub <- subset(clinical, select = c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis",
                                            "clinical_stage", "tumor_grade", "residual_disease_largest_nodule"))

# Merging dataframes by patient as clinical data is made per patient 
core_survival <- merge(clinical_sub, core_survival, by.x = "bcr_patient_barcode", by.y = "X_PATIENT")

## Adding CIBERSORT data
# Stripping last character from sample names
ciber_perc$sample <- sub(".$", "", row.names(ciber_perc))

# Merging data frames
core_survival <- merge(core_survival, ciber_perc, by = "sample")

# Exchange [Not available] to real NAs
core_survival[core_survival == "[Not Available]"] <- NA

data <- core_survival


### [3] Preparing for survival analysis
# Hiloquant function to divide data into low and high
hiloquant<-function(x,pr) {
  qu<-as.vector(quantile(x,pr,na.rm=TRUE))
  a<-which(x>qu)
  b<-which(x<=qu)
  x[a]<-2
  x[b]<-1
  return(x)
}

# Preparing header for saving files
header <- paste("Species", "HR", "loCI", "hiCI", "P", sep = "\t")

# Encoding the data for cox analysis
data_encoded <- data
data_encoded$clinical_stage <- ifelse(data_encoded$clinical_stage == "Stage IC", 1,
                                      ifelse(data_encoded$clinical_stage %in% c("Stage IIA", "Stage IIB", "Stage IIC"), 2,
                                             ifelse(data_encoded$clinical_stage %in% c("Stage IIIa", "Stage IIIB", "Stage IIIC"), 3,
                                                    ifelse(data_encoded$clinical_stage == "Stage IV", 4, NA))))

data_encoded$tumor_grade <- ifelse(data_encoded$tumor_grade == "G1", 11,
                                   ifelse(data_encoded$tumor_grade == "G2", 12,
                                          ifelse(data_encoded$tumor_grade == "G3", 13,
                                                 ifelse(data_encoded$tumor_grade == "G4", 14, 
                                                        ifelse(data_encoded$tumor_grade == "GB", 15, NA)))))

data_encoded$residual_disease_largest_nodule <- ifelse(data_encoded$residual_disease_largest_nodule == "No Macroscopic disease", 20,
                                                       ifelse(data_encoded$residual_disease_largest_nodule == "1-10 mm", 21,
                                                              ifelse(data_encoded$residual_disease_largest_nodule == "11-20 mm", 22,
                                                                     ifelse(data_encoded$residual_disease_largest_nodule == ">20 mm", 23, NA))))



### [4] Survival analysis with immune cells
## Preparing a text file for output
write.table(header, file= file.path(pathout, "/Survival_Immune_cells.txt"), append = TRUE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

## Looping through each cell type
cells <- colnames(ciber_perc)
cells <- cells[cells != "sample"]

for (i in cells) {
  OS_time <- data_encoded$OS.time/30.42
  EVENT <- data_encoded$OS
  ABU <- data_encoded[ , i]
  ex <- hiloquant(as.numeric(ABU), 0.5)
  cov_age <- as.numeric(data_encoded$age_at_initial_pathologic_diagnosis)
  cov_stage <- factor(c(data_encoded$clinical_stage))
  cov_grade <- factor(c(data_encoded$tumor_grade))
  cov_residual <- factor(c(data_encoded$residual_disease_largest_nodule))
  
  # Creating index where clinical parameters are never NA
  index <- which(!is.na(cov_age) & !is.na(cov_stage) & !is.na(cov_grade) & !is.na(cov_residual))
  
  # Running cox with multiple covariates
  cox_res <- coxph(Surv(as.numeric(OS_time[index]), as.numeric(EVENT[index])) ~ ex[index] +
                     cov_age[index] + cov_residual[index] + cov_stage[index] + cov_grade[index])
  
  # Getting summaries
  summary_cox_res <- summary(cox_res)
  coefficient_table <- summary_cox_res$coefficients
  coefficient_of_interest <- coefficient_table["ex[index]", ]
  conf_int <- confint(cox_res)["ex[index]", ]
  
  # Extracting p-value, HR, and CI
  P <- summary_cox_res$coefficients["ex[index]", "Pr(>|z|)"]
  HR <- exp(coefficient_of_interest["coef"])
  loHR <- exp(conf_int[["2.5 %"]])
  upHR <- exp(conf_int[["97.5 %"]])
  loCI <- sprintf("%.2f",round(loHR,2))
  hiCI <- sprintf("%.2f",round(upHR,2))
  
  # Creating text and saving as txt
  text <-paste(i, HR, loCI, hiCI, P, sep="\t")
  #write.table(text, file= file.path(pathout, "/Survival_Immune_cells.txt"),append = TRUE, quote = FALSE, sep = "\t", 
   #           eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  
  # Getting parameters for plotting
  mx<-140
  fit <- survfit(Surv(OS_time[index], EVENT[index]) ~ ex[index])
  logrank <- survdiff(Surv(OS_time[index], EVENT[index]) ~ ex[index])
  n1<-as.vector(logrank[[1]][1])
  m1<-as.vector(logrank[[1]][2])
  
  # Plotting Kaplan Meier Curve
  png(paste0(pathout, "/Plots/KM_Immune_cells/", i, ".png"), width = 20, height = 15, units = "cm", res = 800)
  par(mar = c(9, 5, 4, 2))
  plot(fit, main = paste("Kaplan-Meier Curve for", i),
       xlab = "Time (months)", ylab = "Survival Probability",
       col = c("dodgerblue2", "firebrick"), lwd = 2)
  legend("topright", legend = c("Low", "High"), col = c("dodgerblue2", "firebrick"), lwd = 2)
  HR_text <- sprintf("HR = %.2f", HR)
  P_text <- sprintf("p = %.4f", P)
  text(143, 0.7, HR_text, adj = c(0.5, 0.5), col = "black", cex = 1)
  text(143, 0.6, P_text, adj = c(0.5, 0.5), col = "black", cex = 1)
  OSN1<-length(which(ex[index] == 1))
  OSN2<-length(which(ex[index] == 2))
  timeseq<-seq(20,mx,by=20)
  text(-mx*0.1,-0.3,"No. at risk",xpd=TRUE,col="black",cex=0.8,pos=1)
  text(-mx*0.1,-0.35,"Low abund:",xpd=TRUE,col="dodgerblue2",cex=0.8,pos=1)
  text(-mx*0.1,-0.4,"High abund:",xpd=TRUE,col="firebrick",cex=0.8,pos=1)
  text(0,-0.35,n1,xpd=TRUE,col="dodgerblue2",cex=0.8,pos=1)
  text(0,-0.4,m1,xpd=TRUE,col="firebrick",cex=0.8,pos=1)
  for (j in timeseq) {
    nr1<-n1-length(which(OSN1<j))
    nr2<-m1-length(which(OSN2<j))
    text(j,-0.35,nr1,xpd=TRUE,col="dodgerblue2",cex=0.8,pos=1)
    text(j,-0.4,nr2,xpd=TRUE,col="firebrick",cex=0.8,pos=1)
  }
  dev.off()
}

## Reading in produced table
cox_immune <- read.table(file.path(pathout, "/Survival_Immune_cells.txt"), header = TRUE, sep = "\t")

# Filtering for significant immune cells
cox_immune_fil <- cox_immune[cox_immune$P < 0.01, ]
sig_immune_cells <- cox_immune_fil$Species

cat("Immune cells with significant impact on survival: ", sig_immune_cells)


#### [5] Survival Analysis for individual species
### Preparing a text file for output
write.table(header, file= file.path(pathout, "/Survival_Species.txt"), append = TRUE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

### Looping through each species
species <- core_tax$Names
species <- c("Ernakulamia tanakae", "Verticillium dahliae")

for (i in species) {
  OS_time <- data_encoded$OS.time/30.42
  EVENT <- data_encoded$OS
  ABU <- data_encoded[ , i]
  ex <- hiloquant(as.numeric(ABU), 0.5)
  cov_age <- as.numeric(data_encoded$age_at_initial_pathologic_diagnosis)
  cov_stage <- factor(c(data_encoded$clinical_stage))
  cov_grade <- factor(c(data_encoded$tumor_grade))
  cov_residual <- factor(c(data_encoded$residual_disease_largest_nodule))
  cov_M2 <- as.numeric(c(data_encoded$Macrophages.M2))
  cov_M1 <- as.numeric(c(data_encoded$Macrophages.M1))
  
  # Creating index where clinical parameters are never NA
  index <- which(!is.na(cov_age) & !is.na(cov_stage) & !is.na(cov_grade) & !is.na(cov_residual))
  
  # Running cox with multiple covariates
  cox_res <- coxph(Surv(as.numeric(OS_time[index]), as.numeric(EVENT[index])) ~ ex[index] +
                     cov_age[index] + cov_residual[index] + cov_stage[index] + cov_grade[index] +
                     cov_M2[index] + cov_M1[index])
  
  # Getting summaries
  summary_cox_res <- summary(cox_res)
  coefficient_table <- summary_cox_res$coefficients
  coefficient_of_interest <- coefficient_table["ex[index]", ]
  conf_int <- confint(cox_res)["ex[index]", ]
  
  # Extracting p-value, HR, and CI
  P <- summary_cox_res$coefficients["ex[index]", "Pr(>|z|)"]
  HR <- exp(coefficient_of_interest["coef"])
  loHR <- exp(conf_int[["2.5 %"]])
  upHR <- exp(conf_int[["97.5 %"]])
  loCI <- sprintf("%.2f",round(loHR,2))
  hiCI <- sprintf("%.2f",round(upHR,2))
  
  # Creating text and saving as txt
#  text <-paste(i, HR, loCI, hiCI, P, sep="\t")
#  write.table(text, file= file.path(pathout, "/Survival_Species.txt"),append = TRUE, quote = FALSE, sep = "\t", 
#              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  
  # Getting parameters for plotting
  fit <- survfit(Surv(OS_time[index], EVENT[index]) ~ ex[index])
  logrank <- survdiff(Surv(OS_time[index], EVENT[index]) ~ ex[index])
  n1<-as.vector(logrank[[1]][1])
  m1<-as.vector(logrank[[1]][2])
  
  # Plotting Kaplan Meier Curve
  png(paste0(pathout, "/Plots/KM_Species/", i, ".png"), width = 17, height = 11.5, units = "cm", res = 800)
  par(mar = c(4, 4, 2, 2))
  plot(fit, main = paste(i),
       xlab = "Time (months)", ylab = "Survival Probability",
       col = c("dodgerblue2", "firebrick"), lwd = 2)
  legend("topright", legend = c("Low", "High"), col = c("dodgerblue2", "firebrick"), lwd = 2)
  HR_text <- sprintf("HR = %.2f", HR)
  P_text <- sprintf("p = %.3f", P)
  n1_text <- paste0("No. low = ", n1)
  m1_text <- paste0("No. high = ", m1)
  text(143, 0.75, HR_text, adj = c(0.5, 0.5), col = "black", cex = 1)
  text(143, 0.7, P_text, adj = c(0.5, 0.5), col = "black", cex = 1)
  text(140, 0.6, n1_text, col="dodgerblue2",cex=1)
  text(140, 0.55, m1_text,col="firebrick",cex=1)
  
  dev.off()
}

## Reading in produced table
cox_species <- read.table(file.path(pathout, "/Survival_Species.txt"), header = TRUE, sep = "\t")

# Filtering for significant species
cox_species_fil <- cox_species[cox_species$P < 0.05, ]

sig_species <- cox_species_fil$Species
cat("Species with significant impact on survival: ", sig_species)



#### [6] Survival analysis for bacteria associated with the significant immune cells
## Preparing correlation data for further analysis
correlation <- read.csv(paste0(pathout, "/Correlation_Matrix_fil.csv"), row.names = 1)

# Get species negatively correlating with interesting immune cells
M2_pos <- row.names(correlation[correlation$M2.Macrophages > 0, ])
M2_neg <- row.names(correlation[correlation$M2.Macrophages < 0, ])
M1_pos <- row.names(correlation[correlation$M1.Macrophages > 0, ])
M1_neg <- row.names(correlation[correlation$M1.Macrophages < 0, ])
MemB_pos <- row.names(correlation[correlation$Memory.B.cells > 0, ])
MemB_neg <- row.names(correlation[correlation$Memory.B.cells < 0, ])
DC_pos <- row.names(correlation[correlation$Activated.dendritic.cells > 0, ])
DC_neg <- row.names(correlation[correlation$Activated.dendritic.cells < 0, ])

# Summing species correlating with immune cells
subsets <- c("M2_pos", "M2_neg", "MemB_pos", "MemB_neg") #, "DC_pos", "DC_neg")

for (i in subsets) {
  j <- get(i)
  species_to_keep <- j
  sub_phylo <- subset_taxa(core, rownames(otu_table(core)) %in% j)
  assign(paste0(i, "_phylo"), sub_phylo)
  
  sub_otu <- as.data.frame(otu_table(sub_phylo))
  otu_sums <- colSums(sub_otu)
  sub_otu <- rbind(sub_otu, otu_sums)
  name <- paste0(i, "_sum")
  rownames(sub_otu)[nrow(sub_otu)] <- name
  assign(paste0(i, "_otu"), sub_otu)
  
  sub_tax <- as.data.frame(tax_table(sub_phylo))
  assign(paste0(i, "_tax"), sub_tax)
}


sum_data <- rbind(M2_pos_otu[nrow(M2_pos_otu), ], 
                  M2_neg_otu[nrow(M2_neg_otu), ],
                  #DC_pos_otu[nrow(DC_pos_otu), ],
                  #DC_neg_otu[nrow(DC_neg_otu), ],
                  MemB_pos_otu[nrow(MemB_pos_otu), ],
                  MemB_neg_otu[nrow(MemB_neg_otu), ])


sum_t <- as.data.frame(t(sum_data))

# Creating ratios of positive and negative correlated species
sum_t$M2_ratio <- sum_t$M2_pos_sum / sum_t$M2_neg_sum
sum_t$DC_ratio <- sum_t$DC_pos_sum / sum_t$DC_neg_sum
sum_t$MemB_ratio <- sum_t$MemB_pos_sum / sum_t$MemB_neg_sum

# Log transforming the data
sum_log <- as.data.frame(log(sum_t + 1))

# Stripping last character to have an identical column for merge
sum_log$sample <- sub(".$", "", row.names(sum_log))

# Combining sum data with other data
data_sums <- merge(data_encoded, sum_log, by = "sample")


### Preparing a text file for output
write.table(header, file= file.path(pathout, "/Survival_Corr_Groups.txt"), append = TRUE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

### Looping through the sums
groups <- paste0(subsets, "_sum") 
groups <- c(groups, "M2_ratio", "MemB_ratio") #, "DC_ratio")

for (i in groups) {
  OS_time <- data_sums$OS.time/30.42
  EVENT <- data_sums$OS
  ABU <- data_sums[ , i]
  ex <- hiloquant(as.numeric(ABU), 0.5)
  cov_age <- as.numeric(data_sums$age_at_initial_pathologic_diagnosis)
  cov_stage <- factor(c(data_sums$clinical_stage))
  cov_grade <- factor(c(data_sums$tumor_grade))
  cov_residual <- factor(c(data_sums$residual_disease_largest_nodule))
  cov_M2 <- as.numeric(c(data_sums$Macrophages.M2))
  cov_M1 <- as.numeric(c(data_sums$Macrophages.M1))
  
  # Creating index where clinical parameters are never NA
  index <- which(!is.na(cov_age) & !is.na(cov_stage) & !is.na(cov_grade) & !is.na(cov_residual))
  
  # Running cox with multiple covariates
  cox_res <- coxph(Surv(as.numeric(OS_time[index]), as.numeric(EVENT[index])) ~ ex[index] +
                     cov_age[index] + cov_residual[index] + cov_stage[index] + cov_grade[index] +
                     cov_M2[index] + cov_M1[index])
  
  # Getting summaries
  summary_cox_res <- summary(cox_res)
  coefficient_table <- summary_cox_res$coefficients
  coefficient_of_interest <- coefficient_table["ex[index]", ]
  conf_int <- confint(cox_res)["ex[index]", ]
  
  # Extracting p-value, HR, and CI
  P <- summary_cox_res$coefficients["ex[index]", "Pr(>|z|)"]
  HR <- exp(coefficient_of_interest["coef"])
  loHR <- exp(conf_int[["2.5 %"]])
  upHR <- exp(conf_int[["97.5 %"]])
  loCI <- sprintf("%.2f",round(loHR,2))
  hiCI <- sprintf("%.2f",round(upHR,2))
  
  # Creating text and saving as txt
  text <-paste(i, HR, loCI, hiCI, P, sep="\t")
  write.table(text, file= file.path(pathout, "/Survival_Corr_Groups.txt"),append = TRUE, quote = FALSE, sep = "\t", 
              eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  
  # Plotting Kaplan Meier
  fit <- survfit(Surv(OS_time[index], EVENT[index]) ~ ex[index])
  png(paste0(pathout, "/Plots/KM_Corr_Groups/", i, ".png"), width = 20, height = 15, units = "cm", res = 800)
  plot(fit, main = paste("Kaplan-Meier Curve for", i),
       xlab = "Time (months)", ylab = "Survival Probability",
       col = c("blue", "red"), lwd = 2)
  legend("topright", legend = c("Low", "High"), col = c("blue", "red"), lwd = 2)
  HR_text <- sprintf("HR = %.2f", HR)
  P_text <- sprintf("p = %.4f", P)
  text(140, 0.8, HR_text, adj = c(0.5, 0.5), col = "black", cex = 1)
  text(140, 0.7, P_text, adj = c(0.5, 0.5), col = "black", cex = 1)
  dev.off()
}

## Reading in produced table
cox_groups <- read.table(file.path(pathout, "/Survival_Corr_Groups.txt"), header = TRUE, sep = "\t")

# Filtering for significant species
cox_groups_fil <- cox_groups[cox_groups$P < 0.05, ]

sig_groups <- cox_groups_fil$Species

cat("Groups with significant impact on survival: ", sig_groups)










