library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(neuroCombat)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/Methylation Data - Anil + Lasya/adni.rds')
adni_mri_1 <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/T1_ASHSvols_MTTHK_20241108.csv')
adni_mri_2 <- read.csv('/Users/lasyasreepada/Projects/hippoage/data/ADNI/MRI/T1_ASHSvols_MTTHK_20250221.csv')
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
spareba <- istaging[istaging$Study=='ADNI',]
rm(istaging)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Transform columns
coldata$DateDrawn <- as.Date(coldata$DateDrawn)
dx$EXAMDATE <- as.Date(dx$EXAMDATE)

coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age)

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

# Merge individual spreadsheets
adni_mri <- bind_rows(adni_mri_1, adni_mri_2) %>%
  select(intersect(names(adni_mri_1), names(adni_mri_2)))

# Format metadata
meta <- read.csv('~/Projects/hippoage/data/ADNI/MRI/MR_Image_Acquisition/MRIMETA_All_16Feb2025.csv')
columns <- c('VISCODE', 'VISCODE2')
meta <- meta %>%
  select(-c(X,PHASE,PTID)) %>%
  dplyr::rename(MRIDATE = EXAMDATE) %>%
  mutate(across(all_of(columns), 
                ~str_replace_all(., 
                                 c("scmri" = "bl", 
                                   "blmri" = "bl", 
                                   "sc" = "bl"))))
# Format date columns
adni_mri$MRIDATE <- as.Date(adni_mri$MRIDATE)
meta$MRIDATE <- as.Date(meta$MRIDATE)

# Merge with metadata
adni_mri <- adni_mri %>%
  left_join(meta,by=c('RID','MRIDATE')) %>%
  rowwise() %>%
  mutate(M_Hippo_VOL_ASHST1 = sum(c(M_AHippo_VOL_ASHST1_3T,M_PHippo_VOL_ASHST1_3T)))

# Merge with DNA
adni <- adni %>%
  left_join(adni_mri, by=c('RID'), relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DateDrawn, MRIDATE) %>% as.numeric('years')),
         Age_MRI = abs(interval(PTDOB, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# Select CU only
cu <- adni[adni$DIAGNOSIS==1,]

# Remove missing features
cu <- cu %>%
  filter(!is.na(M_Hippo_VOL_ASHST1))

# Select cross sectional dataset
cu <- cu %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# Harmonize by field strength
dat <- cu %>%
  select(c(ICV, M_Hippo_VOL_ASHST1))
rownames(dat) <- cu$RID
dat <- t(dat)

batch <- cu$FIELD_STRENGTH

age <- cu$Age_MRI
sex <- cu$Sex

mod <- model.matrix(~age+sex)

combat.harmonized <- neuroCombat(dat=dat, batch=batch, mod=mod, parametric=FALSE)
dat.combat <- as.data.frame(t(combat.harmonized$dat.combat))
dat.combat$RID <- as.integer(rownames(dat.combat))

dat.combat <- dat.combat %>%
  dplyr::rename(ICV.combat = ICV,
                M_Hippo_VOL_ASHST1.combat = M_Hippo_VOL_ASHST1)
rownames(dat.combat) <- NULL

# Merge back into main dataframe
cu <- cu %>%
  left_join(dat.combat,by='RID')

# Adjust hippocampal volume for Age & ICV
model <- lm(M_Hippo_VOL_ASHST1.combat ~ ICV.combat, data = cu)
coefficients <- coef(model)

cu$M_Hippo_VOL_ASHST1.adjusted <- cu$M_Hippo_VOL_ASHST1.combat - (coefficients["ICV.combat"] * cu$ICV.combat)

# Select samples
samples <- cu$Sample
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
cu <- left_join(cu,cell_proportions,by="Sample")
cu <- left_join(cu,leukocytes,by='Sample')

# PREPARE FINAL SE
betas_cu <- assays(se)[[1]]
coldata_cu <- as.data.frame(cu)
rownames(coldata_cu) <- coldata_cu$Sample

# Create SE
se <- SummarizedExperiment(assays = betas_cu, colData=coldata_cu)

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
library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")

cg_ok <- !rownames(se) %in% xloci
se <- se[cg_ok,]

# SAVE
saveRDS(se,'~/Projects/hippoage/data/ADNI/Methylation/20250305_HV_SE_ADNI.rds')
