library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ggVennDiagram)
library(SummarizedExperiment)

# READ DATA
aibl <- readRDS("/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl.rds")

# QCS
qcs <- metadata(aibl)
qcs_df <- do.call(rbind, lapply(qcs, as.data.frame))

# DENSITY
ggplot(qcs_df, aes(x=frac_dt)) + geom_density() + theme_pubclean() +
  labs(title="Signal Detection Frequency in AIBL DNA Methylation Samples") +
  xlim(0.9, 1)

cat(sprintf("Mean Signal Detection Frequency: %f\n", mean(qcs_df$frac_dt)))
cat(sprintf("Median Signal Detection Frequency: %f\n", median(qcs_df$frac_dt)))

# LIMIT CPGS

# Read list of common cpgs
cpgs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/cpgs_v1_v2_shared.csv')
cpgs <- cpgs$CPG

# Get betas
betas <- assays(aibl)[[1]]
betas <- betas[cpgs,]

# Get pheno
pheno <- as.data.frame(colData(aibl))

# CREATE SUMMARIZED EXPERIMENT
aibl <- SummarizedExperiment(assays = betas, colData=pheno, metadata = qcs)

# SAVE
saveRDS(aibl, "/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl_new.rds")





