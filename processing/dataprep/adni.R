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
cdr <- read.csv('~/Projects/CoMeT/data/ADNI/cognition/CDR_14Apr2025.csv') # CDR
amyloid <- read.csv('~/Downloads/ADSP_ADNI_PET-Scalar_Dec2023/ADSP_PHC_PET_Amyloid_Simple_Dec2023.csv') # Amyloid Summary
phc <- read.csv('~/Projects/CoMeT/data/ADNI/cognition/ADSP_PHC_COGN_Dec2023.csv') # Cog
ashs_1 <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/T1_ASHSvols_MTTHK_20241108.csv') # ASHS  
ashs_2 <- read.csv('/Users/lasyasreepada/Projects/hippoage/data/ADNI/MRI/T1_ASHSvols_MTTHK_20250221.csv') # ASHS

# READ ISTAGING
# Setup SQLite driver
sqlite.driver <- dbDriver("SQLite")

# Establish database connection
file <- "~/Projects/hippoage/data/istaging.db"
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
  select(c(PTID,Date,APOE_Genotype,APOE4_Alleles,Education_Years,Race,Ethnicity,SPARE_BA)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(SPAREDATE = Date)

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

# Merge individual spreadsheets
ashs <- bind_rows(ashs_1, ashs_2) %>%
  select(intersect(names(ashs_1), names(ashs_2)))

# Format metadata
meta <- read.csv('~/Projects/hippoage/data/ADNI/MRI/MR_Image_Acquisition/MRIMETA_All_16Feb2025.csv')
meta <- meta %>%
  select(c(RID,FIELD_STRENGTH,EXAMDATE)) %>%
  dplyr::rename(MRIDATE = EXAMDATE)

# Format date columns
ashs$MRIDATE <- as.Date(ashs$MRIDATE)
meta$MRIDATE <- as.Date(meta$MRIDATE)

# Merge ASHS with metadata
ashs <- ashs %>%
  left_join(meta,by=c('RID','MRIDATE')) %>%
  rowwise() %>%
  mutate(M_Hippo_VOL_ASHST1 = sum(c(M_AHippo_VOL_ASHST1_3T,M_PHippo_VOL_ASHST1_3T)))

# Merge with DNA
adni <- adni %>%
  left_join(ashs, by=c('RID'), relationship = 'many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DateDrawn, MRIDATE) %>% as.numeric('years')),
         Age_MRI = abs(interval(PTDOB, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# MERGE SPARE

# Merge
adni <- adni %>%
  left_join(spareba, by=c('PTID'), relationship = 'many-to-many') %>%
  mutate(MRItoSPARE = abs(interval(MRIDATE, SPAREDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(MRItoSPARE, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# Select CU only
cu <- adni[adni$DIAGNOSIS==1,]

# Remove missing features
cu <- cu %>%
  filter(if_all(c(M_Hippo_VOL_ASHST1), ~ !is.na(.x)))

# Select cross sectional dataset
cu <- cu %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# HARMONIZATION

# Select features
dat <- cu %>%
  select(c(ICV, M_Hippo_VOL_ASHST1, M_BA35_MeanTHK_MSTTHKMT_ASHST1_3T))
dat <- as.data.frame(dat)
rownames(dat) <- cu$RID
dat <- t(dat)
batch <- cu$FIELD_STRENGTH
age <- cu$Age_MRI
sex <- cu$Sex

# Compose model matrix
mod <- model.matrix(~age+sex)

# Run ComBat
combat.harmonized <- neuroCombat(dat=dat, batch=batch, mod=mod, parametric=FALSE)
dat.combat <- as.data.frame(t(combat.harmonized$dat.combat))
dat.combat$RID <- as.integer(rownames(dat.combat))

# Rename variables
dat.combat <- dat.combat %>%
  dplyr::rename(ICV.combat = ICV,
                M_Hippo_VOL_ASHST1.combat = M_Hippo_VOL_ASHST1,
                M_BA35_MeanTHK_MSTTHKMT_ASHST1.combat = M_BA35_MeanTHK_MSTTHKMT_ASHST1_3T)
rownames(dat.combat) <- NULL

# Merge back into main dataframe
cu <- cu %>%
  left_join(dat.combat,by='RID')

# ADJUST VOLUMES

# Version 1: ICV only
model <- lm(M_Hippo_VOL_ASHST1.combat ~ ICV.combat, data = cu)
coefficients <- coef(model)
cu$M_Hippo_VOL_ASHST1_Adj1 <- cu$M_Hippo_VOL_ASHST1.combat - (coefficients["ICV.combat"] * (cu$ICV.combat - mean(cu$ICV.combat)))

# Version 2: Age and ICV
model2 <- lm(M_Hippo_VOL_ASHST1.combat ~ Age_MRI + ICV.combat, data = cu)
coefficients2 <- coef(model2)
cu$M_Hippo_VOL_ASHST1_Adj2 <- cu$M_Hippo_VOL_ASHST1.combat - (coefficients2["ICV.combat"] * (cu$ICV.combat - mean(cu$ICV.combat)) + coefficients2["Age_MRI"] * (cu$Age_MRI - mean(cu$Age_MRI)))

# ADJUST THICKNESS

# Version 1: Age only
model <- lm(M_BA35_MeanTHK_MSTTHKMT_ASHST1.combat ~ Age_MRI, data = cu)
coefficients <- coef(model)
cu$M_BA35_MeanTHK_MSTTHKMT_ASHST1_Adj <- cu$M_BA35_MeanTHK_MSTTHKMT_ASHST1.combat - (coefficients["Age_MRI"] * (cu$Age_MRI - mean(cu$Age_MRI)))

# SE PREPARATION
# Select samples
samples <- cu$Sample
se <- se[,samples]

# CELL TYPE DECONVOLUTION
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
cu <- left_join(cu,cell_proportions,by="Sample")

# PREPARE SE
betas_cu <- assays(se)[[1]]
coldata_cu <- as.data.frame(cu)
rownames(coldata_cu) <- coldata_cu$Sample

# Create SE
se <- SummarizedExperiment(assays = betas_cu, colData=coldata_cu)

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
cleanMatrixForClusterSE <- function(se, f_row = 0, f_col = 0.25) {
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

# Get betas
betas <-  assays(se)[[1]]
betas <- as.data.frame(betas)

# Convert row names to a column
betas <- betas %>% 
  tibble::rownames_to_column("CpG")

# Get pheno
coldata <- as.data.frame(colData(se))

# SAVE
library(arrow)
write_parquet(betas,'/Users/lasyasreepada/Projects/hippoage/data/ADNI/20250427_adni_betas.parquet')
write.csv(coldata,'/Users/lasyasreepada/Projects/hippoage/data/ADNI/20250427_adni_pheno.csv')

saveRDS(se,'~/Projects/hippoage/data/ADNI/20250612_ADNI_SE.rds')
