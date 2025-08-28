library(glmnet)
library(caret)
library(sesame)
library(pracma)
library(SummarizedExperiment)
library(doParallel)
library(arrow)

sesameDataCache()

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/Methylation/20250228_HV_SE.rds')

betas <- assays(se)[[1]]
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Function to impute missing CpGs
imputeRowMean <- function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  return(mtx)
}

# Impute missing CpGs
betas <- imputeRowMean(betas)

# Create dataframe
data <- as.data.frame(t(betas))
data$Sex <- as.integer(coldata$Sex == "Female")
data$Leuk <- as.integer(coldata$Leuk >= 0.5)
data$HV <- coldata$M_Hippo_VOL_ASHST1.adjusted

# Write
write_parquet(data, "~/Projects/hippoage/data/20250228_ADNI_EN.parquet")

# Load data (assuming it's in a data frame called 'data')
se <- readRDS('~/Projects/hippoage/data/ADNI/Methylation/20250227_SPAREBA_SE.rds')
se$SPAREBAG <- se$SPARE_BA.SVM.RBF. - se$Age_MRI

betas <- assays(se)[[1]]
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Function to impute missing CpGs
imputeRowMean <- function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  return(mtx)
}

# Impute missing CpGs
betas <- imputeRowMean(betas)

# Create dataframe
data <- as.data.frame(t(betas))
data$Sex <- as.integer(coldata$Sex == "Female")
data$Leuk <- as.integer(coldata$Leuk >= 0.5)
data$SPAREBAG <- coldata$SPAREBAG

# Write
write_parquet(data, "~/Projects/hippoage/data/20250228_SPAREBA_EN.parquet")
