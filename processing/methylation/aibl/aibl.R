library(dplyr)
library(readxl)
library(sesame)
library(BiocParallel)
library(SummarizedExperiment)

sesameDataCache("EPIC.address")
sesameDataCache("KYCG.EPIC.Mask.20220123")

# READ DATA
idat_dir <- "/Users/lasyasreepada/Projects/hippoage/data/AIBL/GSE153712_RAW/"
betas <- openSesame(idat_dir, BPPARAM = BiocParallel::MulticoreParam(2))
qcs <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats, funs="detection",BPPARAM = BiocParallel::MulticoreParam(2))
pvals <- openSesame(idat_dir, func = pOOBAH, return.pval=TRUE,BPPARAM = BiocParallel::MulticoreParam(2))

# SAVE
save(betas, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/betas.RData')
save(qcs, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/qcs.RData')
save(pvals, file = '/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/pvals.RData')

# MERGE WITH ANNOTATION
pheno <- read_excel('/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Corey_age_project_DNAmAge_data_cognition_imaging_20240916.xlsx')

# Rename some columns
pheno <- pheno %>% 
  dplyr::rename(AIBL_ID = AIBL_ID...1,
                Sample_Title = MWAS_ID) %>%
  select (-c(...31, ...32,AIBL_ID...33))

# Extract the sample names
samples <- colnames(betas)
part1 <- sapply(strsplit(samples, "_"), `[`, 2)
part2 <- sapply(strsplit(samples, "_"), `[`, 3)
samples_abridged <- paste(part1, part2, sep = "_")

# Create a dataframe
samples <- as.data.frame(samples)
names(samples) <- 'Sample'
samples$Sample_Title <- samples_abridged

# Merge
pheno <- left_join(samples,pheno,by='Sample_Title')

# CREATE SUMMARIZED EXPERIMENT
aibl <- SummarizedExperiment(assays = betas, colData=pheno, metadata = qcs)

# SAVE
saveRDS(aibl, "/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl.rds")


