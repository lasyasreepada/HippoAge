library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(RSQLite)
library(readxl)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/Processed/INDD_SE_20250201.rds') # DNA 
clocks <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/Processed/indd_clocks_20250201.csv') # Clocks 
dx <- read_excel('~/Projects/hippoage/data/INDD/NC.xlsx')
cdr <- read_excel('~/Projects/hippoage/data/INDD/CDR.xlsx')
apoe <- read_excel('~/Projects/hippoage/data/INDD/APOE.xlsx')
amyloid <- read.csv('~/Projects/hippoage/data/INDD/abc_multimodal_amystatus_2025March.csv')
ashs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_MRI3T_ASHST1T2_measurements_20241220.csv')

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

# Select PENN
spare <- istaging[istaging$Study=='PENN',]
rm(istaging)
spare <- spare %>%
  select(c(PTID,Date,SPARE_BA)) %>%
  dplyr::rename(INDDID = PTID,
                SPAREBA = SPARE_BA,
                SPAREDATE = Date)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Transform columns
coldata$DNAdate <- as.Date(coldata$DNAdate)

coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age)

# MERGE CLOCKS
clocks <- clocks %>%
  select(-c(X,Sample_Well,Sample_Plate,DNAdate,Version,INDDID,DOB,Sex,Age))
coldata <- coldata %>%
  left_join(clocks, by = c("Sample"))

# MERGE DX
dx <- dx %>%
  dplyr::rename(DXDATE = VisitDate)
dx$DXDATE <- as.Date(dx$DXDATE)
dx$INDDID <- as.character(dx$INDDID)

# Merge
abc <- coldata %>%
  left_join(dx, by = c("INDDID"),relationship = 'many-to-many') %>%
  mutate(DNAtoDX = abs(interval(DNAdate, DXDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(INDDID, DNAdate) %>% # Group by ID and Date1
  slice_min(DNAtoDX, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE CDR

# Format
cdr <- cdr %>%
  dplyr::rename(CDRDATE = TestDate)
cdr$CDRDATE <- as.Date(cdr$CDRDATE)
cdr$INDDID <- as.character(cdr$INDDID)

# Merge
abc <- abc %>%
  left_join(cdr, by = c("INDDID"),relationship = 'many-to-many') %>%
  mutate(DNAtoCDR = abs(interval(DNAdate, CDRDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(INDDID, DNAdate) %>% # Group by ID and Date1
  slice_min(DNAtoCDR, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE APOE
apoe <- apoe %>%
  dplyr::rename(APOE_Genotype = APOE)
apoe$INDDID <- as.character(apoe$INDDID)
apoe$APOE4_Alleles <- str_count(apoe$APOE_Genotype, "4")

abc <- abc %>%
  left_join(apoe, by = c("INDDID"))

# MERGE AMYLOID
amyloid <- amyloid %>%
  dplyr::rename(AMYDATE = AmyDate)
amyloid$INDDID <- as.character(amyloid$INDDID)
amyloid$AMYDATE <- as.Date(amyloid$AMYDATE,format = "%m/%d/%Y")
amyloid$Amy_Status <- as.integer(amyloid$Amy_Status == "Positive")

abc <- abc %>%
  left_join(amyloid, by = c("INDDID"),relationship = 'many-to-many') %>%
  mutate(DNAtoAMY = abs(interval(DNAdate, AMYDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(INDDID, DNAdate) %>% # Group by ID and Date1
  slice_min(DNAtoAMY, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE MRI
# Format date columns
ashs <- ashs %>%
  dplyr::rename(MRIDATE = SCANDATE,
                INDDID = INDD)
ashs$INDDID <- as.character(ashs$INDDID)
ashs$MRIDATE <- as.Date(as.character(ashs$MRIDATE), format = "%Y%m%d")

# Merge with DNA
abc <- abc %>%
  left_join(ashs, by=c('INDDID'), relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DNAdate, MRIDATE) %>% as.numeric('years')),
         Age_MRI = abs(interval(DOB, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(INDDID, DNAdate) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# MERGE SPARE

# Format
spare$SPAREDATE <- as.Date(spare$SPAREDATE)
spare$INDDID <- as.character(spare$INDDID)
spare$INDDID <- gsub("_", ".", spare$INDDID)

# Merge
abc <- abc %>%
  left_join(spare, by=c('INDDID'), relationship = 'many-to-many') %>%
  mutate(MRItoSPARE = abs(interval(MRIDATE, SPAREDATE) %>% as.numeric('years'))) %>%
  group_by(INDDID, DNAdate) %>%
  slice_min(MRItoSPARE, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# Select CU only
abc <- abc[abc$NORMCOG==1,]

# Remove missing features
abc <- abc %>%
  filter(!is.na(SPAREBA))

# Select cross sectional dataset
abc <- abc %>%
  group_by(INDDID) %>%
  slice_min(MRItoSPARE, with_ties = FALSE) %>%
  ungroup()

# Select samples
samples <- abc$Sample
se <- se[,samples]

# CELL TYPE DECONVOLUTION
# install_github("zhou-lab/knowYourCG")
betas <- assays(se)[[1]]

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
abc <- left_join(abc,cell_proportions,by="Sample")

# PREPARE FINAL SE
betas <- assays(se)[[1]]
coldata <- as.data.frame(abc)
rownames(coldata) <- coldata$Sample

# Create SE
se <- SummarizedExperiment(assays = betas, colData=coldata)

# Limit to CpGs in common between EPIC v1 and EPIC v2
common_cpgs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/cpgs_v1_v2_shared.csv')
common_cpgs <- common_cpgs$CPG
se <- se[common_cpgs,]

# Filter cross-reactive probes
library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")

cg_ok <- !rownames(se) %in% xloci
se <- se[cg_ok,]

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
saveRDS(se,'~/Projects/hippoage/data/INDD/20250312_ABC_SE.rds')
