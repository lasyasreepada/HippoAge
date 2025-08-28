library(dplyr)
library(sesame)
library(BiocParallel)
library(SummarizedExperiment)

sesameDataCache("EPICv2.address")
sesameDataCache("KYCG.EPICv2.Mask.20230314")
sesameDataCache("EPIC.address")
sesameDataCache("KYCG.EPIC.Mask.20220123")

# DATA IOs
data_dir <- "/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation"
idat_dirs <- list.dirs(path=data_dir, full.names = TRUE, recursive = TRUE)
idat_subfolders <- idat_dirs[grepl("/20", idat_dirs)]

# TEST ONE DIRECTORY
total_items <- length(idat_subfolders)

# RUN SESAME (EPIC v2)

# First batch to initialize matrices
idat_dir <- idat_subfolders[1]
cat(sprintf("Processing batch %s\n", basename(idat_dir)))
betas_v2 <- openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(2))
qcs_v2 <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats, funs="detection",BPPARAM = BiocParallel::MulticoreParam(2))
pvals_v2 <- openSesame(idat_dir, func = pOOBAH, return.pval=TRUE,BPPARAM = BiocParallel::MulticoreParam(2))
cat(sprintf("Finished batch %d of %d\n", 1, total_items))

for (i in 2:85){
  # Process each subfolder
  idat_dir <- idat_subfolders[i]
  cat(sprintf("Processing batch %s\n", basename(idat_dir)))
  betas_i <- openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(2))
  qcs_i <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats, funs="detection",BPPARAM = BiocParallel::MulticoreParam(2))
  pvals_i <- openSesame(idat_dir, func = pOOBAH, return.pval=TRUE,BPPARAM = BiocParallel::MulticoreParam(2))
  cat(sprintf("Finished batch %d of %d\n", i, total_items))
  
  # Combine outputs
  betas_v2 <- cbind(betas_v2, betas_i)
  qcs_v2 <- c(qcs_v2, qcs_i)
  pvals_v2 <- cbind(pvals_v2, pvals_i)
  
}

# WRITE BATCH 1
save(betas_v2, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/betas_EPICv2.RData')
save(qcs_v2, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/qcs_EPICv2.RData')
save(pvals_v2, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/pvals_EPICv2.RData')

# RUN SESAME (EPIC v1)

# First batch to initialize matrices
idat_dir <- idat_subfolders[86]
cat(sprintf("Processing batch %s\n", basename(idat_dir)))
betas_v1 <- openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(2))
qcs_v1 <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats, funs="detection",BPPARAM = BiocParallel::MulticoreParam(2))
pvals_v1 <- openSesame(idat_dir, func = pOOBAH, return.pval=TRUE,BPPARAM = BiocParallel::MulticoreParam(2))
cat(sprintf("Finished batch %d of %d\n", 86, total_items))

for (i in 87:total_items){
  # Process each subfolder
  idat_dir <- idat_subfolders[i]
  cat(sprintf("Processing batch %s\n", basename(idat_dir)))
  betas_i <- openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(2))
  qcs_i <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats, funs="detection",BPPARAM = BiocParallel::MulticoreParam(2))
  pvals_i <- openSesame(idat_dir, func = pOOBAH, return.pval=TRUE,BPPARAM = BiocParallel::MulticoreParam(2))
  cat(sprintf("Finished batch %d of %d\n", i, total_items))
  
  # Combine outputs
  betas_v1 <- cbind(betas_v1, betas_i)
  qcs_v1 <- c(qcs_v1, qcs_i)
  pvals_v1 <- cbind(pvals_v1, pvals_i)
  
}

# WRITE BATCH 2
save(betas_v1, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/betas_EPICv1.RData')
save(qcs_v1, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/qcs_EPICv1.RData')
save(pvals_v1, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/pvals_EPICv1.RData')

# MERGE WITH ANNOTATION
pheno <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/MethIDs_8May2023.csv')

pheno <- pheno %>% 
  dplyr::rename(ID = Saple_Nae,
         Sample = MethylomeID)

samples_v1 <- as.data.frame(colnames(betas_v1))
samples_v2 <- as.data.frame(colnames(betas_v2))

names(samples_v1) <- 'Sample'
names(samples_v2) <- 'Sample'

pheno_v1 <- left_join(samples_v1,pheno,by='Sample')
pheno_v2 <- left_join(samples_v2,pheno,by='Sample')

# CREATE SUMMARIZED EXPERIMENT
abc_v1 <- SummarizedExperiment(assays = betas_v1, colData=pheno_v1, metadata = qcs_v1)
abc_v2 <- SummarizedExperiment(assays = betas_v2, colData=pheno_v2, metadata = qcs_v2)
  
# SAVE
saveRDS(abc_v1, "/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/abc_v1.rds")
saveRDS(abc_v2, "/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_Methylation_Processed/abc_v2.rds")

# Some random stuff
# https://www.bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/inferences.html
# identical(colnames(assay(se)), rownames(colData(se)))

