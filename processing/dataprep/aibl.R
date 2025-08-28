library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(RSQLite)
library(neuroCombat)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Methylation/Processed/aibl.rds')
ashs <- read.csv('/Users/lasyasreepada/Projects/hippoage/data/AIBL/MRI/T1_ASHSvols_MTTHK_20250319_TEMP.csv')
meta <- read.csv('/Users/lasyasreepada/Projects/hippoage/data/AIBL/MRI/AIBL_FieldStrength_20250331.csv')

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
spareba <- istaging[istaging$Study=='AIBL',]
rm(istaging)

# Select required columns
spareba <- spareba %>%
  select(c(PTID,Date,SPARE_BA)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(SPAREDATE = Date)

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

# Transform columns
coldata <- coldata %>%
  dplyr::rename(Age_DNA = Age,
                RID = AIBL_ID)

# Get DOB
date_str <- paste0(as.character(coldata$Demographic.YearMonthOfBirth), "01")
coldata$DOB <- as.Date(date_str, format = "%Y%m%d")
ashs$MRIDATE <- as.Date(ashs$MRIDATE)

# MERGE MRI
ashs <- ashs %>%
  rowwise() %>%
  mutate(M_Hippo_VOL_ASHST1 = sum(c(M_AHippo_VOL_ASHST1_3T,M_PHippo_VOL_ASHST1_3T)))

# Format metadata
meta <- meta %>%
  select(c(Subject.ID,Study.Date,Imaging.Protocol)) %>%
  dplyr::rename(RID = Subject.ID,
                MRIDATE = Study.Date,
                Field_Strength = Imaging.Protocol)

meta$Field_Strength <- as.numeric(gsub(".*=(\\d+\\.\\d+)", "\\1", meta$Field_Strength))

# Format date columns
ashs$MRIDATE <- as.Date(ashs$MRIDATE)
meta$MRIDATE <- as.Date(meta$MRIDATE,  format = "%m/%d/%Y")

# Merge ASHS with metadata
ashs <- ashs %>%
  left_join(meta,by=c('RID','MRIDATE'))

# Merge with DNA
aibl <- coldata %>%
  left_join(ashs, by=c('RID'), relationship = 'many-to-many') %>%
  mutate(Age_MRI = abs(interval(DOB, MRIDATE) %>% as.numeric('years')),
         DNAtoMRI = abs(Age_DNA - Age_MRI)) %>%
  group_by(RID, Age_DNA) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# MERGE SPARE
spareba$RID <- as.integer(spareba$PTID)
spareba <- spareba %>%
  select(-c(PTID))

# Merge
aibl <- aibl %>%
  left_join(spareba, by=c('RID'), relationship = 'many-to-many') %>%
  mutate(MRItoSPARE = abs(interval(MRIDATE, SPAREDATE) %>% as.numeric('years'))) %>%
  group_by(RID, Age_DNA) %>%
  slice_min(MRItoSPARE, with_ties = FALSE) %>%
  ungroup()

# PREPARE DATASET

# Select aibl only
aibl <- aibl[aibl$Diagnosis=="HC",]

# Remove missing features
aibl <- aibl %>%
  filter(if_all(c(M_Hippo_VOL_ASHST1,Field_Strength), ~ !is.na(.x)))

# Select cross sectional dataset
aibl <- aibl %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# HARMONIZATION

# Select features
dat <- aibl %>%
  select(c(ICV, M_Hippo_VOL_ASHST1))
dat <- as.data.frame(dat)
rownames(dat) <- aibl$RID
dat <- t(dat)
batch <- aibl$Field_Strength
age <- aibl$Age_MRI
sex <- aibl$Sex

# Compose model matrix
mod <- model.matrix(~age+sex)

# Run ComBat
combat.harmonized <- neuroCombat(dat=dat, batch=batch, mod=mod, parametric=FALSE)
dat.combat <- as.data.frame(t(combat.harmonized$dat.combat))
dat.combat$RID <- as.integer(rownames(dat.combat))

# Rename variables
dat.combat <- dat.combat %>%
  dplyr::rename(ICV.combat = ICV,
                M_Hippo_VOL_ASHST1.combat = M_Hippo_VOL_ASHST1)
rownames(dat.combat) <- NULL

# Merge back into main dataframe
aibl <- aibl %>%
  left_join(dat.combat,by='RID')

# ADJUST VOLUMES

# Version 1: ICV only
model <- lm(M_Hippo_VOL_ASHST1.combat ~ ICV.combat, data = aibl)
coefficients <- coef(model)
aibl$M_Hippo_VOL_ASHST1_Adj1 <- aibl$M_Hippo_VOL_ASHST1.combat - (coefficients["ICV.combat"] * (aibl$ICV.combat - mean(aibl$ICV.combat)))

# Version 2: Age and ICV
model <- lm(M_Hippo_VOL_ASHST1.combat ~ Age_MRI + ICV.combat, data = aibl)
coefficients <- coef(model)
aibl$M_Hippo_VOL_ASHST1_Adj2 <- aibl$M_Hippo_VOL_ASHST1.combat - (coefficients["ICV.combat"] * (aibl$ICV.combat - mean(aibl$ICV.combat)) + coefficients["Age_MRI"] * (aibl$Age_MRI - mean(aibl$Age_MRI)))

# For Pipeline
aibl$M_BA35_MeanTHK_MSTTHKMT_ASHST1_Adj <- aibl$M_BA35_MeanTHK_MSTTHKMT_ASHST1_3T
  
# Select samples
samples <- aibl$Sample
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
aibl <- left_join(aibl,cell_proportions,by="Sample")

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

# Get Pheno
coldata <- as.data.frame(colData(se))

# SAVE
library(arrow)
write_parquet(betas,'/Users/lasyasreepada/Projects/hippoage/data/AIBL/20250427_aibl_betas.parquet')
write.csv(coldata,'/Users/lasyasreepada/Projects/hippoage/data/AIBL/20250427_aibl_pheno.csv')
