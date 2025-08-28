library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggVennDiagram)
library(SummarizedExperiment)

# READ DATA
indd <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/Processed/INDD_SE_20250201.rds')

# QCS
qcs <-metadata(indd)
qcs_df <- do.call(rbind, lapply(qcs, as.data.frame))

# DENSITY
ggplot(qcs_df, aes(x=frac_dt)) + geom_density() + theme_pubclean() +
  labs(title="Signal Detection Frequency in INDD DNA Methylation Samples")

mean(qcs_df$frac_dt)
median(qcs_df$frac_dt)

# OLD
# abc_v1 <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/abc_v1.rds')
# abc_v2 <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/abc_v2.rds')
# qcs_v1 <- metadata(abc_v1)
# qcs_v2 <- metadata(abc_v2)
# qcs_v1_df <- do.call(rbind, lapply(qcs_v1, as.data.frame))
# qcs_v2_df <- do.call(rbind, lapply(qcs_v2, as.data.frame))
# p1 <- ggplot(qcs_v1_df, aes(x=frac_dt)) + geom_density() + theme_pubclean() +
#   labs(title="Signal Detection Frequency in ABC DNA Methylation Samples - EPICv1") +
#   xlim(0.2, 1)
# 
# p2 <- ggplot(qcs_v2_df, aes(x=frac_dt)) + geom_density() + theme_pubclean() +
#   labs(title="Signal Detection Frequency in ABC DNA Methylation Samples - EPICv2") +
#   xlim(0.2, 1)
# 
# p1 + p2 + plot_layout(ncol = 1)
# mean(qcs_v1_df$frac_dt)
# mean(qcs_v2_df$frac_dt)
# median(qcs_v1_df$frac_dt)
# median(qcs_v2_df$frac_dt)

# VENN

# BETAS
betas_v1 <- assays(abc_v1)[[1]]
betas_v2 <- assays(abc_v2)[[1]]

# CPG names
cpgs_v1 <- rownames(betas_v1)
cpgs_v2 <- rownames(betas_v2)

# Get rid of tail
cpgs_v2 <- sapply(strsplit(cpgs_v2, "_"), `[`, 1)

# Venn list
venn_data <- list(v1 = cpgs_v1, v2 = cpgs_v2)

# Plot the Venn diagram
ggVennDiagram(venn_data) +
  scale_fill_gradient(
    low = "lightblue",  # Color for low values
    high = "darkblue",  # Color for high values
  ) +
  ggtitle('Venn Diagram of CpGs across EPICv1 and EPICv2')
