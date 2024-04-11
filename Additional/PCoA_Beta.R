#### Beta diversity
pathout_OV <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/OV_Data_nt"
pathout_PAAD <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja/PAAD_Data_nt"
pathout <- "/data/projects/2020/OvarianCancerHH/Thesis_Katja"


library(phyloseq)


micro_OV <- readRDS(paste0(pathout_OV, "/microorganisms_decontam.rds"))
micro_PAAD <- readRDS(paste0(pathout_PAAD, "/microorganisms_decontam.rds")) 

micro <- merge_phyloseq(micro_OV, micro_PAAD)

### Ordination using PCoA and bray distance
ord_PCOA <- ordinate(micro, "PCoA")

options(repr.plot.width = 3, repr.plot.height = 3)
plot_ordination(micro, ord_PCOA, type = "samples", color = "project.x") +
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "line")) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("TCGA-PAAD" = "tomato", "TCGA-OV" = "slateblue"), name = "TCGA project") +
  xlab("PCoA [32.5%]") + ylab("PCoA [11.6%]")

ggsave(paste0(pathout, "/Plots/Beta_PCoA.png"), height = 8, width = 12, unit = "cm")


