library(SummarizedExperiment)
library(dplyr)

# READ DATA
se <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/20250218_HV_SE.rds')
cg_scan <- read.csv('~/Projects/hippoage/data/ADNI/EWAS/20250218_HV_CGs.csv')

# Slice to include cgs from scan
cg_probes <- cg_scan$Probe_ID
se <- se[cg_probes,]

# Get betas
betas <- assays(se)[[1]]

# Get coldata
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Operations

# Impute missing CpGs
imputeRowMean <- function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  return(mtx)
}

betas <- imputeRowMean(betas)
mx <- t(betas)
mx <- as.data.frame(mx)

# # Compute variance for each CpG site
# cpg_variances <- apply(betas, 2, var)
# 
# # Filter CpGs with variance above a threshold
# threshold <- quantile(cpg_variances, 0.9)  # Keep top 10% most variable
# betas <- betas[, cpg_variances >= threshold]

# Residualize hippocampal volume for Age & ICV
lm_hv <- lm(M_Hippo_VOL_ASHST1_3T ~ Age + ICV, data = coldata)

# Add Features
mx$HV_res <- residuals(lm_hv)
mx$Sex <- coldata$Sex

write.csv(mx,'~/Projects/hippoage/data/ADNI/20250218_EN_Data.csv')



