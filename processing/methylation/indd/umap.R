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

# indd <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/Processed/INDD_SE_20250201.rds')

indd <- readRDS('~/Projects/hippoage/data/INDD/20250312_ABC_SE.rds')

# CLEAN
indd <- cleanMatrixForClusterSE(indd)

# BETAS
mx <- assays(indd)[[1]]

# PARTICIPANTS
pheno <- as.data.frame(colData(indd))
pheno <- pheno[!is.na(pheno$ID),]
ids <- pheno$ID
samples <- pheno[pheno$ID %in% ids,'Sample']

# Slice Betas Matrix
# mx <- mx[ ,colnames(mx) %in% samples]

# MATRIX OPERATIONS
mx <- imputeRowMean(mx)
mx <- t(mx)

# UMAP
library(umap)
umap <- umap(mx)
umap_layout <- as.data.frame(umap$layout)
colnames(umap_layout) <- c("UMAP1", "UMAP2")
pheno <- cbind(as_tibble(pheno),umap_layout)

write.csv(pheno, '~/Projects/hippoage/data/INDD/umap.csv')


p1 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Sex)) + geom_point() +
  labs(title='Sex') # Check on two outlying probes

# ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Age)) + geom_point() +
#   labs(title='Age') # Check on two outlying probes

p2 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Race)) + geom_point() +
  labs(title='Race') 

# ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Ethnicity)) + geom_point() +
#   labs(title='Ethnicity') # Check on two outlying probes

# ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Deceased)) + geom_point() +
#   labs(title='Deceased') 

p3 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Version)) + geom_point() +
  labs(title='EPIC Version')

p4 <- ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Sample_Plate)) + geom_point() +
  labs(title='Sample Plate')

# ggplot(pheno,aes(x=UMAP1,y=UMAP2,color=Sample_Well)) + geom_point(show.legend = FALSE) +
#   labs(title='Sample_Well')

grid.arrange(p1, p2, p3, p4, nrow = 2)









