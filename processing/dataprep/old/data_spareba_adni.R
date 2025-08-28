library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(RSQLite)
library(maxprobes)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/Methylation Data - Anil + Lasya/adni.rds')
clocks <- read.csv('~/Projects/CoMeT/data/ADNI/methylation/clocks.csv')
dx <- read.csv('~/Projects/hippoage/data/ADNI/DXSUM_21Jan2025.csv')

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

# Select ADNI
adni_mri <- istaging[istaging$Study=='ADNI',]
rm(istaging)
adni_mri <- adni_mri %>%
  select(c(PTID,Phase,Date,SPARE_BA.SVM.RBF.,Age)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(Age_MRI = Age,
                MRIDATE = Date)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Transform columns
coldata$DateDrawn <- as.Date(coldata$DateDrawn)
dx$EXAMDATE <- as.Date(dx$EXAMDATE)

coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age)

# MERGE CLOCKS
clocks <- clocks %>%
  select(-c(RID,Phase,Edate,DateDrawn,PlateNumber,Array,Slide,PTDOB,Sex,Age,X))
coldata <- coldata %>%
  left_join(clocks, by = c("Sample"))

# MERGE DX
dx <- dx %>%
  select(c(RID,PTID,VISCODE2,EXAMDATE,DIAGNOSIS))

adni <- coldata %>%
  left_join(dx, by = c("RID"),relationship = 'many-to-many') %>%
  mutate(DNAtoDX = abs(interval(DateDrawn, EXAMDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(RID, DateDrawn) %>% # Group by ID and Date1
  slice_min(DNAtoDX, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE MRI

# Format date columns
adni_mri$MRIDATE <- as.Date(adni_mri$MRIDATE)

# Merge with DNA
adni <- adni %>%
  left_join(adni_mri, by=c('PTID'), relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DateDrawn, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# SELECT CU
adni <- adni[adni$DIAGNOSIS==1,]

# Remove missing features
adni <- adni %>%
  filter(!is.na(SPARE_BA.SVM.RBF.))

# Select cross sectional dataset
adni <- adni %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# Select samples
samples <- adni$Sample
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
adni <- left_join(adni,cell_proportions,by="Sample")
adni <- left_join(adni,leukocytes,by='Sample')

# PREPARE FINAL SE
betas_adni <- assays(se)[[1]]
coldata_adni <- as.data.frame(adni)
rownames(coldata_adni) <- coldata_adni$Sample

# Create SE
se <- SummarizedExperiment(assays = betas_adni, colData=coldata_adni)

# Limit to CpGs in common between EPIC v1 and EPIC v2
common_cpgs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/cpgs_v1_v2_shared.csv')
common_cpgs <- common_cpgs$CPG
se <- se[common_cpgs,]

# MISSING VALUE HANDLING

# Function to remove rows and features with too much missingness
cleanMatrixForClusterSE <- function(se, f_row = 0.25, f_col = 0.25) {
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

# Filter
se <- cleanMatrixForClusterSE(se)

# Categorical level checks
colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
cg_ok <- (checkLevels(assay(se), colData(se)$Sex))

# Filter
se <- se[cg_ok,]

# Filter cross-reactive probes
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")

cg_ok <- !rownames(se) %in% xloci
se <- se[cg_ok,]

# SAVE
saveRDS(se,'~/Projects/hippoage/data/ADNI/EWAS/20250305_SPAREBA_SE_ADNI.rds')
