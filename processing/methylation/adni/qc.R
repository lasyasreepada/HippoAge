library(ggplot2)
library(Rtsne)
library(sesame)
library(stats)
library(factoextra)
library(dplyr)
library(matrixStats)
library(ExperimentHub)
library(SummarizedExperiment)

setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read SE
dnam <- readRDS("data/ADNI/methylation/methylation.rds")

# Extract coldata as a dataframe
coldata <- as.data.frame(colData(dnam))

# Read adni dataframe
adni <- read.csv("data/ADNI/MRI_DNA_COG_X_Grouped.csv")

# Additional information in coldata

# Extract betas
betas <- assays(adni)[[1]]

# Extract QC
qcs <- metadata(adni)
adniQC <- do.call(rbind, lapply(qcs, as.data.frame))

# Density Plots
ggplot(adniQC, aes(x=frac_dt)) + geom_density() + theme_pubclean() +
  labs(title="Signal Detection Frequency in DNA Methylation Samples")

# Additional QCs
idat_dir_adni <- "/Users/sreepada/Library/CloudStorage/Box-Box/DataRepository/ADNI/Methylation/ADNI_iDAT_files/"
qcs_intensity <- openSesame(idat_dir_adni, prep="", func=sesameQC_calcStats, funs="intensity")
adniQC_intensity <- do.call(rbind, lapply(qcs_intensity, as.data.frame))

