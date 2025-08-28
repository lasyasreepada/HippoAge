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
se <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/20250305_SPAREBA_SE_ADNI.rds')
se$SPAREBA <- se$SPARE_BA.SVM.RBF.

coldata <- as.data.frame(colData(se))

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
# SPARE BA
res <- DML(se,~Age_MRI+Sex+PC1_Blood+PC2_Blood+SPAREBA, BPPARAM=BiocParallel::MulticoreParam(6))
# res <- DML(se,~Age_MRI+Sex+Leuk+SPAREBA, BPPARAM=BiocParallel::MulticoreParam(6))

# Test
smry <- summaryExtractTest(res)
smry$FDR <- p.adjust(smry$Pval_SPAREBA, method="fdr")
smry %>% arrange(Pval_SPAREBA) %>% dplyr::select(Probe_ID,Est_SPAREBA,Pval_SPAREBA,FDR)

# SAVE 
saveRDS(smry, '~/Projects/hippoage/data/ADNI/EWAS/20250305_SPAREBA_EWAS2.rds')




