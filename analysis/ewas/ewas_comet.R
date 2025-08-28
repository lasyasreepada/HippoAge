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
se <- readRDS('~/Projects/hippoage/data/ADNI/20250319_COMET_SE.rds')
se <- se[,se$Case==1]
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# RUN EWAS

# CoMeT
res <- DML(se,~PlateNumber+CD4Tnv+Baso+CD4Tmem+Bmem+Treg+CD8Tmem+CD8Tnv+Eos+NK+Neu+Mono+CorticalLOAD_CoMeT, BPPARAM=BiocParallel::MulticoreParam(6))
res <- DML(se,~Age_MRI+PlateNumber+CD4Tnv+Baso+CD4Tmem+Bmem+Treg+CD8Tmem+CD8Tnv+Eos+NK+Neu+Mono+CorticalLOAD_CoMeT, BPPARAM=BiocParallel::MulticoreParam(6))

# Test
smry <- summaryExtractTest(res)
smry$FDR <- p.adjust(smry$Pval_CorticalLOAD_CoMeT, method="fdr")
smry %>% arrange(Pval_CorticalLOAD_CoMeT) %>% dplyr::select(Probe_ID,Est_CorticalLOAD_CoMeT,Pval_CorticalLOAD_CoMeT,FDR) %>%
  filter(Pval_CorticalLOAD_CoMeT < 1e-5)

# SAVE 
saveRDS(smry, '~/Projects/hippoage/data/ADNI/EWAS/20250319_COMET_SMRY_NOAGE.rds')
saveRDS(smry, '~/Projects/hippoage/data/ADNI/EWAS/20250319_COMET_SMRY_AGE.rds')
