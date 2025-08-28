library(SummarizedExperiment)
library(dplyr)
library(Rtsne)
library(ggplot2)
library(gridExtra)
library(scales)

#remove rows and features with too much missingness
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

#impute missing cgs
imputeRowMean <- function(mtx) {
    k <- which(is.na(mtx), arr.ind=TRUE)
    mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
    return(mtx)
}

tsneSE <- function(se, perplexity=30, seed=1) {
    library(Rtsne)
    set.seed(seed)
    se = cleanMatrixForClusterSE(se)
    mx = imputeRowMean(assay(se))
    ## samples = colnames(mx)
    tsne = Rtsne(t(mx), dims=2, perplexity=perplexity)
    df = as.data.frame(tsne$Y)
    colnames(df) = c("tSNE1", "tSNE2")
    df$sample = colnames(mx)
    cbind(df, as_tibble(colData(se))) #[samples,]))
}

pcaSE <- function(se) {
    se = cleanMatrixForClusterSE(se, f_row=0.7, f_col=0.7)
    mx = imputeRowMean(assays(se)[[1]])
    samples = colnames(mx)
    pca = prcomp(t(mx))
    cbind(pca$x, as_tibble(colData(se)[samples,]))
    ## cbind(as_tibble(colData(se)), pca$x[samples,])
    ## cbind(as_tibble(colData(se)), pca$x[colData(se)$IDAT,]) # or this
}

umapSE <- function(se, seed=1) {
    library(umap)
    set.seed(seed)
    se = cleanMatrixForClusterSE(se, f_row=0.7, f_col=0.7)
    mx = imputeRowMean(assay(se))
    samples = colnames(mx)
    umapResult <- data.frame(umap(t(mx))$layout)
    colnames(umapResult) <- c("UMAP1", "UMAP2")
    cbind(umapResult, as_tibble(colData(se)[samples,]))
}

# Read SE and ADNI
se <- readRDS('~/Projects/hippoage/data/ADNI/20250612_ADNI_SE.rds')
se <- cleanMatrixForClusterSE(se)

# Extract betas as matrix
mx <- assays(se)[[1]]
coldata <- as.data.frame(colData(se))

# Select participants of interest
# adni <- read.csv('~/Projects/CoMeT/data/ADNI/MRI_DNA_COG_X_Grouped.csv')
# adni <- adni[!is.na(adni$AGE_SAMPLE_nearest),]
# rids <- adni$RID
# samples <- coldata[coldata$RID %in% rids,'Sample']
# mx <- mx[ ,colnames(mx) %in% samples]

# Matrix operations
mx <- imputeRowMean(mx)
mxt <- t(mx)

# UMAP
library(umap)
umapResult <- umap(mxt)
umapResult2 <- as.data.frame(umapResult$layout)
colnames(umapResult2) <- c("UMAP1", "UMAP2")
coldata <- cbind(umapResult2, as_tibble(coldata))

# Rename some columns
coldata <- coldata %>%
  rename(Diagnosis = DIAGNOSIS,
         Plate = PlateNumber)

# Factor some columns
coldata$Diagnosis <- as.factor(coldata$Diagnosis)
coldata$Plate <- as.factor(coldata$Plate)

# Write coldata for easy access later:
write.csv(coldata,'~/Projects/CoMeT/data/ADNI/methylation/umap.csv')

# Plots
p1 <- ggplot(coldata,aes(x=UMAP1,y=UMAP2,color=Sex)) + geom_point() +
  labs(title='Sex')
p3 <- ggplot(coldata,aes(x=UMAP1,y=UMAP2,color=Array)) + geom_point() +
  labs(title='Array')
p4 <- ggplot(coldata,aes(x=UMAP1,y=UMAP2,color=Plate)) + geom_point() + guides(colour = guide_legend(nrow = 10)) +
  labs(title='Plate')
p2 <- ggplot(coldata,aes(x=UMAP1,y=UMAP2,color=Race)) + geom_point() + scale_color_hue(direction = -1) +
  labs(title='Race')

grid.arrange(p1, p2, p3, p4, nrow = 2)

# TSNE
library(Rtsne)
set.seed(1)
tsne <- Rtsne(mxt, dims=2, perplexity=30,verbose = TRUE, max_iter = 100)

# MDS plot
library(minfi)

mdsPlot(
  mx
)
# df_pca <- pcaSE(se)
# df_tsne <- tsneSE(se)
# df_umap <- umapSE(se)
# ggplot(df,aes(x=PC1,y=PC2,color=Pathology_Group)) + geom_point()













