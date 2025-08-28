library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(RSQLite)
library(maxprobes)
library(neuroCombat)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/Methylation Data - Anil + Lasya/adni.rds') # DNA
clocks <- read.csv('~/Projects/CoMeT/data/ADNI/methylation/clocks.csv') # Clocks
dx <- read.csv('~/Projects/hippoage/data/ADNI/DXSUM_21Jan2025.csv') # Diagnosis
cdr <- read.csv('~/Downloads/Neuropsychological/CDR_09Mar2025.csv') # CDR
amyloid <- read.csv('~/Downloads/ADSP_ADNI_PET-Scalar_Dec2023/ADSP_PHC_PET_Amyloid_Simple_Dec2023.csv') # Amyloid Summary
phc <- read.csv('~/Projects/CoMeT/data/ADNI/cognition/ADSP_PHC_COGN_Dec2023.csv') # Cog
fs <- read.csv('/Users/lasyasreepada/Projects/CoMeT/data/ADNI/imaging/UCSFFSX51_11_08_19_14Apr2025.csv')

# READ ISTAGING
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
spareba <- istaging[istaging$Study=='ADNI',]
rm(istaging)

# Select required columns
spareba <- spareba %>%
  select(c(PTID,APOE_Genotype,APOE4_Alleles,Education_Years,Race,Ethnicity)) %>%
  filter(complete.cases(.)) %>%
  group_by(PTID) %>%
  slice_head() %>%
  ungroup()

spareba$Education_Years <- as.integer(spareba$Education_Years)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Transform coldata
coldata$DateDrawn <- as.Date(coldata$DateDrawn)
coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age)

# MERGE CLOCKS
clocks <- clocks %>%
  select(-c(RID,Phase,Edate,DateDrawn,PlateNumber,Array,Slide,PTDOB,Sex,Age,X))
coldata <- coldata %>%
  left_join(clocks, by = c("Sample"))

# MERGE DX

# Format
dx <- dx %>%
  select(c(RID,PTID,VISCODE2,EXAMDATE,DIAGNOSIS)) %>%
  dplyr::rename(DXDATE = EXAMDATE)
dx$DXDATE <- as.Date(dx$DXDATE)

# Merge
adni <- coldata %>%
  left_join(dx, by = c("RID"),relationship = 'many-to-many') %>%
  mutate(DNAtoDX = abs(interval(DateDrawn, DXDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(RID, DateDrawn) %>% # Group by ID and Date1
  slice_min(DNAtoDX, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE CDR

# Format
cdr <- cdr %>%
  select(c(RID,VISDATE,CDRSB))
cdr$VISDATE <- as.Date(cdr$VISDATE)

# Merge
adni <- adni %>%
  left_join(cdr, by = c("RID"),relationship = 'many-to-many') %>%
  mutate(DNAtoCDR = abs(interval(DateDrawn, VISDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(RID, DateDrawn) %>% # Group by ID and Date1
  slice_min(DNAtoCDR, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE PHC

# Format
phc <- phc %>%
  select(c(RID,EXAMDATE,PHC_Education,PHC_MEM,PHC_EXF,PHC_LAN,PHC_VSP)) %>%
  dplyr::rename(PHCDATE = EXAMDATE)
phc$PHCDATE <- as.Date(phc$PHCDATE)

# Merge 
adni <- adni %>%
  left_join(phc, by = c("RID"),relationship = 'many-to-many') %>%
  mutate(DNAtoPHC = abs(interval(DateDrawn, PHCDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(RID, DateDrawn) %>% # Group by ID and Date1
  slice_min(DNAtoPHC, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE AMYLOID

# Format
amyloid <- amyloid %>% 
  select(c(RID,SCANDATE,PHC_CENTILOIDS,PHC_AMYLOID_STATUS)) %>%
  dplyr::rename(AMYDATE = SCANDATE)
amyloid$AMYDATE <- as.Date(as.character(amyloid$AMYDATE), format = "%Y%m%d")

# Merge 
adni <- adni %>%
  left_join(amyloid, by = c("RID"),relationship = 'many-to-many') %>%
  mutate(DNAtoAMY = abs(interval(DateDrawn, AMYDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(RID, DateDrawn) %>% # Group by ID and Date1
  slice_min(DNAtoAMY, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE MRI
fs <- fs %>%
  select(c(RID,Case,EXAMDATE,CorticalLOAD_CoMeT,CorticalEOAD_CoMeT)) %>%
  dplyr::rename(MRIDATE = EXAMDATE)
fs$MRIDATE <- as.Date(fs$MRIDATE)

# Merge with DNA
adni <- adni %>%
  left_join(fs, by=c('RID'),relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DateDrawn, MRIDATE) %>% as.numeric('years')),
         Age_MRI = abs(interval(PTDOB, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# MERGE SPARE

# Merge
adni <- adni %>%
  left_join(spareba, by=c('PTID'))

# PREPARE DATASET

# Remove missing features
adni_temp <- adni %>%
  filter(if_all(c(CorticalLOAD_CoMeT), ~ !is.na(.x)))

# Select cross sectional dataset
adni_temp <- adni_temp %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# Final Data
adni_comet <- adni_temp[adni_temp$Case==1,]

# SE 

# Select samples
samples <- adni_comet$Sample
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
adni_comet <- left_join(adni_comet,cell_proportions,by="Sample")

# PREPARE FINAL SE
betas <- assays(se)[[1]]
coldata <- as.data.frame(adni_comet)
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

# MISSING VALUE HANDLING

# Function to remove rows and features with too much missingness
cleanMatrixForClusterSE <- function(se, f_row = 0, f_col = 0.5) {
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

# SAVE
saveRDS(se,'~/Projects/hippoage/data/ADNI/20250319_COMET_SE.rds')
