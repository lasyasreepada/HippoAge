library(dplyr)
library(SummarizedExperiment)

# Read MSA data
msa <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/20250528_raw_betas_432.rds')
cpgs <- rownames(msa)
cpgs <- sapply(strsplit(cpgs, "_", fixed = TRUE), `[`, 1)
cpgs <- unique(cpgs)

# Save to file
write.csv(cpgs, '~/Projects/hippoage/data/msa_cpgs.csv',row.names = FALSE)

# CpGs from Clocks
age <- read.csv('~/Projects/hippoage/data/model_coefs_Age_DNA.csv')
spareba <- read.csv('~/Projects/hippoage/data/model_coefs_SPARE_BA.csv')
hv1 <- read.csv('~/Projects/hippoage/data/model_coefs_M_Hippo_VOL_ASHST1_Adj1.csv')

age <- age[(!age$Feature.Names == 'Intercept' & !age$Coefficients==0),]
spareba <- spareba[(!spareba$Feature.Names == 'Intercept' & !spareba$Coefficients==0),]
hv1 <- hv1[(!hv1$Feature.Names == 'Intercept' & !hv1$Coefficients==0),]

# Get overlap
age_overlap <- intersect(cpgs, age$Feature.Names)
spareba_overlap <- intersect(cpgs, spareba$Feature.Names)
hv1_overlap <- intersect(cpgs, hv1$Feature.Names)