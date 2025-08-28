library(dplyr)
library(Rtsne)
library(scales)
library(ggplot2)
library(gridExtra)
library(SummarizedExperiment)

# Remove rows and features with too much missingness
cleanMatrixForClusterSE <- function(se, f_row = 0.5, f_col = 0.5) {
  mtx = assays(se)[[1]]
  cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
              f_row, f_col))
  cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
  namtx = is.na(mtx)
  good_row = rowSums(namtx) <= ncol(mtx) * f_row
  good_col = colSums(namtx) <= nrow(mtx) * f_col
  cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
  return(se[good_row, good_col])
}

# Impute missing CpGs
imputeRowMean <- function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  return(mtx)
}

# READ DATA
aibl <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl_new.rds')

# CLEAN
aibl <- cleanMatrixForClusterSE(aibl)

# BETAS
mx <- assays(aibl)[[1]]

# PHENOTYPES
pheno <- as.data.frame(colData(aibl))

# MATRIX OPERATIONS
mx <- imputeRowMean(mx)
mxt <- t(mx)

# UMAP
library(umap)
umap <- umap(mxt)
umap_layout <- as.data.frame(umap$layout)
colnames(umap_layout) <- c("UMAP1", "UMAP2")
pheno <- cbind(as_tibble(pheno),umap_layout)

# SAVE
write.csv(pheno,'/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/phenotypes_umap.csv')

# PLOTS
p1 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Age)) + geom_point() +
  labs(title='Age')

p2 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Sex)) + geom_point() +
  labs(title='Sex')

p3 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Diagnosis)) + geom_point() +
  labs(title='Diagnosis')

p4 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=as.factor(apoe))) + geom_point() +
  labs(title='APOE4 Alleles')

p5 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Smoking_ALL)) + geom_point() +
  labs(title='Smoking Status')

p6 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=as.factor(Centiloid))) + geom_point() +
  labs(title='Centiloid')

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


