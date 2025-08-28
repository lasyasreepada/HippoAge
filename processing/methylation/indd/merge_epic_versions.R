library(dplyr)
library(SummarizedExperiment)

# READ DATA
abc_v1 <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/Old/abc_v1.rds')
abc_v2 <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/Old/abc_v2.rds')

# BETAS
betas_v1 <- assays(abc_v1)[[1]]
betas_v2 <- assays(abc_v2)[[1]]

# CPGS
cpgs_v1 <- rownames(betas_v1)
cpgs_v2 <- rownames(betas_v2)

# Get rid of tail
cpgs_v2 <- sapply(strsplit(cpgs_v2, "_"), `[`, 1)

# Reset rownames for betas v2
rownames(betas_v2) <- cpgs_v2

# Intersection
cpgs <- intersect(cpgs_v1,cpgs_v2)

# Save CPG list
cpgs <- as.data.frame(cpgs)
names(cpgs) <- 'CPG'
write.csv(cpgs,'/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/cpgs_v1_v2_shared.csv')

# Edit arrays
betas_v1 <- betas_v1[cpgs,]
betas_v2 <- betas_v2[cpgs,]

# MERGE BETAS
betas <- cbind(betas_v1,betas_v2)

# PHENO
pheno_v1 <- as.data.frame(colData(abc_v1))
pheno_v1$Version <- 'EPICv1'

pheno_v2 <- as.data.frame(colData(abc_v2))
pheno_v2$Version <- 'EPICv2'

pheno <- rbind(pheno_v1,pheno_v2)

# QCS 
qc_v1 <- metadata(abc_v1)
qc_v2 <- metadata(abc_v2)
qc <- c(qc_v1,qc_v2)

# CREATE SUMMARIZED EXPERIMENT
abc <- SummarizedExperiment(assays = betas, colData=pheno, metadata = qc)

# SAVE
saveRDS(abc, "/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/indd.rds")

# UPDATE ANNOTATION

# Read 
indd <- readRDS("/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/indd.rds")
annot <- read_excel('/Users/lasyasreepada/Projects/hippoage/data/INDD/MethylationIDs.xlsx')
annot2 <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/MethIDs_8May2023.csv')

library(data.table)

# Convert to data.table
setDT(annot)
setDT(annot2)

annot2 <- annot2 %>%
  dplyr::rename(INDDID= Saple_Nae)

# Update NA values in table1 using table2
pheno <- annot[annot2, on = "INDDID", Methylome1_DataID := ifelse(is.na(Methylome1_DataID), i.MethylomeID, Methylome1_DataID)]
pheno$INDDID <- (sapply(strsplit(pheno$ID, "_"), `[`, 1))

write.csv(ids,'~/Projects/hippoage/data/INDD/IDs.csv')






