library(dplyr)
library(tableone)
library(ggVennDiagram)

aibl_mri <- read.csv('/Users/lasyasreepada/Projects/hippoage/data/AIBL/AIBL_MRI_12_29_2024.csv')
aibl_dna <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl.rds')

# GET COLDATA
pheno <- as.data.frame(colData(aibl_dna))
rownames(pheno) <- NULL

pheno <- pheno %>%
  dplyr::rename(Age_DNA = Age)

# RENAME
aibl_mri <- aibl_mri %>%
  select(-c(Group, Format, Downloaded)) %>%
  dplyr::rename(AIBL_ID = Subject,
                MRIDATE = Acq.Date,
                Age_MRI = Age)

# MERGE
aibl <- pheno %>%
  left_join(aibl_mri, by = c("AIBL_ID"),relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(Age_DNA - Age_MRI)) %>% # Calculate the absolute date difference
  group_by(AIBL_ID, Age_DNA) %>% # Group by ID and Date1
  slice_min(DNAtoMRI, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# Remove missing features
aibl <- aibl %>%
  filter(!is.na(DNAtoMRI))

# SELECT CU
aibl_cu <- aibl[aibl$Diagnosis=='HC',]

# Select cross sectional dataset
aibl_cu <- aibl_cu %>%
  group_by(AIBL_ID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# SAVE 
write.csv(aibl,'/Users/lasyasreepada/Projects/hippoage/data/AIBL/AIBL.csv')



