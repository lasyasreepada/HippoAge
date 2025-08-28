library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(RSQLite)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl.rds')

# Setup SQLite driver
sqlite.driver <- dbDriver("SQLite")

# Establish database connection
file <- "~/Projects/predict4ad/data/istaging.db"
db <- dbConnect(sqlite.driver, dbname=file)

# Read in brain age data
dbListTables(db)
istaging <- dbReadTable(db,"istaging")

# Disconnect SQLite
dbDisconnect(db)

# Select aibl
aibl_mri <- istaging[istaging$Study=='AIBL',]
rm(istaging)
aibl_mri <- aibl_mri %>%
  select(c(PTID,Phase,Date,SPARE_BA.SVM.RBF.,Age)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(Age_MRI = Age,
                MRIDATE = Date)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Transform columns
coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age,
                PTID = AIBL_ID)

aibl_mri$PTID <- as.integer(aibl_mri$PTID)
aibl_mri$MRIDATE <- as.Date(aibl_mri$MRIDATE)

# MERGE MRI

# Merge with DNA
aibl <- coldata %>%
  left_join(aibl_mri, by=c('PTID'), relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(Age_DNA - Age_MRI)) %>%
  group_by(PTID, Age_DNA) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# Select CU only
# aibl <- aibl[aibl$Diagnosis=="HC",]

# Remove missing features
aibl <- aibl %>%
  filter(!is.na(SPARE_BA.SVM.RBF.))

# Select cross sectional dataset
aibl <- aibl %>%
  group_by(PTID) %>%
  slice_min(DNAtoSPARE, with_ties = FALSE) %>%
  ungroup()

# Select samples
samples <- aibl$Sample
se <- se[,samples]

# CELL TYPE DECONVOLUTION
# install_github("zhou-lab/knowYourCG")
betas <- assays(se)[[1]]

# Estimate proportion of leukocytes
leukocytes <- estimateLeukocyte(betas,platform = 'EPIC')
leukocytes <- data.frame(Sample = names(leukocytes), Leuk = as.numeric(leukocytes))
aibl <- left_join(aibl,leukocytes,by='Sample')

# PREPARE FINAL SE
betas_aibl <- assays(se)[[1]]
coldata_aibl <- as.data.frame(aibl)
rownames(coldata_aibl) <- coldata_aibl$Sample

# Create SE
se <- SummarizedExperiment(assays = betas_aibl, colData=coldata_aibl)

# Limit to CpGs in common between EPIC v1 and EPIC v2
common_cpgs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/cpgs_v1_v2_shared.csv')
common_cpgs <- common_cpgs$CPG
se <- se[common_cpgs,]

# # MISSING VALUE HANDLING
# 
# # Function to remove rows and features with too much missingness
# cleanMatrixForClusterSE <- function(se, f_row = 0.25, f_col = 0.25) {
#   mtx = assays(se)[[1]]
#   cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
#               f_row, f_col))
#   cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
#   namtx = is.na(mtx)
#   good_row = rowSums(namtx) <= ncol(mtx) * f_row
#   good_col = colSums(namtx) <= nrow(mtx) * f_col
#   cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
#   return(se[good_row, good_col])
# }
# 
# # Filter
# se <- cleanMatrixForClusterSE(se)
# 
# # Categorical level checks
# colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
# cg_ok <- (checkLevels(assay(se), colData(se)$Sex))
# 
# # Filter
# se <- se[cg_ok,]

# SAVE
saveRDS(se,'~/Projects/hippoage/data/AIBL/AIBL_SPAREBA_SE.rds')
