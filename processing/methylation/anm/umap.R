library(dplyr)
library(ggplot2)
library(ggpubr)
library(GEOquery)
library(SummarizedExperiment)

# READ DATA
gset <- getGEO("GSE144858", GSEMatrix =TRUE, getGPL=FALSE)
gset <- gset[[1]]

methylation_data <- as.data.frame(exprs(gset))
sample_info <- pData(gset)

# CREATE SE
anm <- SummarizedExperiment(assays = list(betas=methylation_data), colData=sample_info)

# SAVE
saveRDS(anm,'~/Projects/hippoage/data/ANM/anm.rds')

# UMAP
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

# CLEAN
anm <- cleanMatrixForClusterSE(anm)

# BETAS
betas <- as.matrix(assays(anm)[[1]])

# PHENOTYPES
pheno <- as.data.frame(colData(anm))

# MATRIX OPERATIONS
betas <- imputeRowMean(betas)
betast <- t(betas)

# UMAP
library(umap)
umap <- umap(betast)
umap_layout <- as.data.frame(umap$layout)
colnames(umap_layout) <- c("UMAP1", "UMAP2")
pheno <- cbind(as_tibble(pheno),umap_layout)

# REFACTOR
pheno$age.ch1 <- as.integer(pheno$age.ch1)

# PLOTS
ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=age.ch1)) + geom_point() +
  labs(title='Age')

ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Sex.ch1)) + geom_point() +
  labs(title='Sex')

ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=disease.state.ch1)) + geom_point() +
  labs(title='Diagnosis')

ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=age.65.exclusion.ch1)) + geom_point() +
  labs(title='Under 65')

# # QC
# qc <- read.csv('~/Downloads/GSE144858_Matrix_Signal_Intensities.csv')
# 
# dt <- grep('_Detection.pval', names(qc), value=TRUE)
# qc <- qc[,dt]
# qc_avg <- qc %>% summarise(across(everything(), mean))
# qc_avg <- as.data.frame(t(qc_avg), col.names='DT')
# qc_avg <- 1 - qc_avg
# rownames(qc_avg) <- NULL
# qc_avg <- t(qc_avg)
# 
# ggplot(qc_avg, aes(x=V1)) + geom_density() + theme_pubclean() +
#   labs(title="Signal Detection Frequency in DNA Methylation Samples") +
#   xlim(0.99,1)
