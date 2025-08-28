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
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/Processed/INDD_SE_20250201.rds')
clocks <- read.csv('~/Projects/hippoage/data/INDD/indd_clocks_20250201.csv')

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

# Select Penn
penn_mri <- istaging[istaging$Study=='PENN',]
rm(istaging)
penn_mri <- penn_mri %>%
  select(c(PTID,Phase,Date,SPARE_BA.SVM.RBF.,Age,Diagnosis_nearest_2.0)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(Age_MRI = Age,
                MRIDATE = Date)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# MERGE CLOCKS
clocks <- clocks %>%
  select(-c(X,Sample_Well,Sample_Plate,DNAdate,Version,INDDID,DOB,Sex,Age))
coldata <- coldata %>%
  left_join(clocks, by = c("Sample"))

# MERGE MRI
# Transform columns
coldata$DNAdate <- as.Date(coldata$DNAdate)

coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age,
                PTID = INDDID)

penn_mri$MRIDATE <- as.Date(penn_mri$MRIDATE)
penn_mri$PTID <- gsub("_", ".", penn_mri$PTID)

# Merge with DNA 
penn <- coldata %>%
  left_join(penn_mri, by=c('PTID'), relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DNAdate, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(PTID, DNAdate) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# Select CU only
# penn <- penn[penn$DIAGNOSIS==1,]

# Remove missing features
penn <- penn %>%
  filter(!is.na(SPARE_BA.SVM.RBF.))

# Select cross sectional dataset
penn <- penn %>%
  group_by(PTID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# Select samples
samples <- penn$Sample
se <- se[,samples]

# CELL TYPE DECONVOLUTION
# install_github("zhou-lab/knowYourCG")
betas <- assays(se)[[1]]

# Estimate proportion of leukocytes
leukocytes <- estimateLeukocyte(betas,platform = 'EPIC')
leukocytes <- data.frame(Sample = names(leukocytes), Leuk = as.numeric(leukocytes))

# Cell type proportions
library(EpiDISH)
library(FlowSorted.Blood.EPIC)
data(cent12CT.m)

# Set reference matrix
ref <- cent12CT.m

# Intersect
common_cpgs <- intersect(rownames(betas), rownames(ref))
betas <- betas[common_cpgs,]
ref <- ref[common_cpgs,]

# Run deconvolution
deconv <- epidish(beta.m = betas, 
                  ref.m = ref, 
                  method = "RPC")

# Extract cell type proportions
cell_proportions <- as.data.frame(deconv$estF)
cell_proportions <- rownames_to_column(cell_proportions,var = "Sample")

# Merge
penn <- left_join(penn,cell_proportions,by="Sample")
penn <- left_join(penn,leukocytes,by='Sample')

# PREPARE FINAL SE
betas_penn <- assays(se)[[1]]
coldata_penn <- as.data.frame(penn)
rownames(coldata_penn) <- coldata_penn$Sample

# Create SE
se <- SummarizedExperiment(assays = betas_penn, colData=coldata_penn)

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
saveRDS(se,'~/Projects/hippoage/data/INDD/PENN_SPAREBA_SE.rds')
