library(dplyr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(ggplot2)
library(ggVennDiagram)
library(tableone)
library(devtools)
library(tibble)
library(doParallel)

sesameDataCache()

# READ DATA
se <- readRDS('~/Projects/hippoage/data/ADNI/Methylation/20250312_HV_SE_ADNI.rds')
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# PCA on Blood Cell Composition
# Load required library
library(stats)

# Perform PCA
blood <- coldata[,c("CD4Tnv","Baso","CD4Tmem","Bmem","Bnv","Treg","CD8Tmem","CD8Tnv","Eos","NK","Neu","Mono")]
pca_result <- prcomp(blood, scale = TRUE)

# Extract the first 2 principal components
pc1 <- pca_result$x[, 1]
pc2 <- pca_result$x[, 2]

se$PC1_Blood <- pc1
se$PC2_Blood <- pc2

# RUN EWAS

# Harmonized Hippocampal Volume
res <- DML(se,~Age_MRI+Sex+Leuk+ICV.combat+M_Hippo_VOL_ASHST1.combat, BPPARAM=BiocParallel::MulticoreParam(6))
res <- DML(se,~Age_MRI+Sex+PC1_Blood+PC2_Blood+ICV.combat+M_Hippo_VOL_ASHST1.combat, BPPARAM=BiocParallel::MulticoreParam(6))
res <- DML(se,~Age_MRI+Sex+CD4Tnv+Baso+CD4Tmem+Bmem+Treg+CD8Tmem+CD8Tnv+Eos+NK+Neu+Mono+ICV.combat+M_Hippo_VOL_ASHST1.combat, BPPARAM=BiocParallel::MulticoreParam(6))

# Test
smry <- summaryExtractTest(res)
smry$FDR <- p.adjust(smry$Pval_M_Hippo_VOL_ASHST1.combat, method="fdr")
smry %>% arrange(Pval_M_Hippo_VOL_ASHST1.combat) %>% dplyr::select(Probe_ID,Est_M_Hippo_VOL_ASHST1.combat,Pval_M_Hippo_VOL_ASHST1.combat,FDR)

# SAVE 
saveRDS(smry, '~/Projects/hippoage/data/ADNI/EWAS/20250308_HV_ADNI_FINAL.rds')

smry <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/20250308_HV_ADNI_FINAL.rds')
smry %>% arrange(Pval_M_Hippo_VOL_ASHST1.combat) %>% dplyr::select(Probe_ID,Est_M_Hippo_VOL_ASHST1.combat,Pval_M_Hippo_VOL_ASHST1.combat,FDR)

