### Coverage plots ###
library(ggplot2)

pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja"
path_cancer <- paste0(pathout, "/PAAD_Data_nt")
path_samtools <- paste0(path_cancer, "/Genome_alignment/1522002_CDS_samtools.coverage")
num_patients <- 179
Species <- "Actinomyces succiniciruminis"

# Reading in data
cds <- read.table(path_samtools,
                     header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Changing column names
colnames(cds) <- c("CDS", "locus", "depth")

# Create new dataframe where values for each CDS can be placed
cds_new <- data.frame(CDS = unique(cds$CDS),
                      Length = 0,
                      kbp = 0,
                      Read_sum = 0,
                      RPK = 0
)

# Extracting length, kbp and sums of reads for each CDS
unique_cds <- unique(cds$CDS)

length <- c()
kbp <- c()
read_sum <- c()

for (i in unique_cds) {
  filter_df <- cds[cds$CDS == i, ]
  length <- c(length, (max(filter_df$locus) - min(filter_df$locus)))
  just_length <- max(filter_df$locus) - min(filter_df$locus)
  kbp <- c(kbp, just_length / 1000)
  read_sum <- c(read_sum, sum(filter_df$depth))
}

# Adding information to the new dataframe
cds_new$Length <- length
cds_new$kbp <- kbp
cds_new$Read_sum <- read_sum
cds_new$RPK <- cds_new$Read_sum / cds_new$kbp
cds_new$RPKM <- (as.numeric(cds_new$Read_sum) * 1e9) / (as.numeric(cds_new$Length) * sum(as.numeric(cds_new$Read_sum)))
cds_new$RPKM_patient <- cds_new$RPKM / num_patients

# Plotting the RPKM values
ggplot(cds_new, aes(x = CDS, y = RPKM_patient)) +
  geom_col(col = "#377eb8") + 
  ylim(0, 90)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank()) +
  ggtitle(Species) +
  xlab("Coding sequence") +
  ylab("RPKM/patient") +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(pathout, "/Plots/Verticillium.png"), height = 6, width = 10, unit = "cm")

# "#377eb8", "#984ea3", "#e41a1c"



