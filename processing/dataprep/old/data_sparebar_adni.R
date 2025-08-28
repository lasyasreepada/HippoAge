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
library(lme4)

sesameDataCache()

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/Methylation Data - Anil + Lasya/adni.rds') # DNA
clocks <- read.csv('~/Projects/CoMeT/data/ADNI/methylation/clocks.csv') # Clocks
dx <- read.csv('~/Projects/hippoage/data/ADNI/DXSUM_21Jan2025.csv') # Diagnosis
cdr <- read.csv('~/Downloads/Neuropsychological/CDR_09Mar2025.csv') # CDR
amyloid <- read.csv('~/Downloads/ADSP_ADNI_PET-Scalar_Dec2023/ADSP_PHC_PET_Amyloid_Simple_Dec2023.csv') # Amyloid Summary
phc <- read.csv('~/Projects/CoMeT/data/ADNI/cognition/ADSP_PHC_COGN_Dec2023.csv') # Cog
ashs_1 <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/T1_ASHSvols_MTTHK_20241108.csv') # ASHS  
ashs_2 <- read.csv('/Users/lasyasreepada/Projects/hippoage/data/ADNI/MRI/T1_ASHSvols_MTTHK_20250221.csv') # ASHS

# ashs_1 <- ashs_1 %>%
#   filter(!is.na(M_Hippo_MeanTHK_MSTTHKMT_ASHST1_3T))
# 
# ashs_2 <- ashs_2[!ashs_2$RID %in% intersect(ashs_1$RID,ashs_2$RID),]
# 
# ashs_2 <- ashs_2 %>%
#   filter(!is.na(M_Hippo_MeanTHK_MSTTHKMT_ASHST1_3T))

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
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# MERGE SPARE

# Merge
adni <- adni %>%
  left_join(spareba, by=c('PTID'), relationship = 'many-to-many') %>%
  mutate(MRItoSPARE = interval(MRIDATE, SPAREDATE) %>% as.numeric('years'))

# DX FILTER
adni <- adni[adni$DIAGNOSIS==1,]

spare <- adni %>%
  group_by(SPAREDATE) %>%
  slice_sample(n=1)

library(data.table)

# Count the number of timepoints
grouped <- as.data.table(spare)[, count := uniqueN(SPARE_BA), by=PTID]

# Extract longitudinal data with > than 2 timepoints overall
adni_spare <- grouped[count >= 2,]

# Filter data for baseline condition (e.g., time = 0 and condition = "A")
baseline_data <- adni_spare %>%
  group_by(RID) %>%
  slice_min(MRItoSPARE, with_ties = FALSE)

# Extract baseline values
baseline_values <- baseline_data %>%
  select(RID, SPARE_BA) %>%
  dplyr::rename(SPARE_BA_BL = SPARE_BA)

# Merge baseline values back into the original dataset
adni_spare <- merge(adni_spare, baseline_values, by = "RID")

# LMEM Model
response <- "SPARE_BA"
# time since baseline + Baseline SPARE_BA, 
lmem_formula = as.formula(paste0(response, " ~ MRItoSPARE + SPARE_BA_BL + (1 + MRItoSPARE|PTID)"))

# Fit LMEM
lmem <- lmer(lmem_formula, data=adni_spare, control=lmerControl(optimizer ="Nelder_Mead",
                                                              calc.derivs = FALSE, optCtrl = list(method = "optimx", starttests = FALSE, kkt = FALSE)))

fixed_effects <- fixef(lmem)
random_effects <- ranef(lmem)$PTID
fixed_slopes <- fixed_effects['MRItoSPARE']
random_slopes <- random_effects['MRItoSPARE']
subject_slopes <- fixed_slopes + random_slopes

names(subject_slopes) <- 'SPARE_BA_slope'
subject_slopes$PTID <- rownames(subject_slopes)
adni_spare <- left_join(adni_spare,subject_slopes,by="PTID")

# PLOTTING
adni_spare$SPARE_BA_pred <- fitted(lmem)

library(ggplot2)
library(ggforce)

# Create plot for selected subjects with the individual regression lines overlaid
myplot <- ggplot(adni_spare) + 
  geom_point(aes(x = MRItoSPARE, y = SPARE_BA), color = "black") +
  geom_smooth(method='lm',aes(x = MRItoSPARE, y = SPARE_BA_pred)) + 
  facet_wrap_paginate(~PTID, nrow=4, ncol = 4, page=1)

# Calculate the number of pages needed
n_pages <- n_pages(myplot)

# Create a multi-page PDF
pdf("plot.pdf", width = 8.5, height = 11)

for (i in 1:n_pages) {
  print(
    myplot + facet_wrap_paginate(~PTID, nrow=4, ncol = 4, page=i)
  )
}

dev.off()

adni_spare <- adni_spare %>%
  group_by(RID) %>%
  slice_min(MRItoSPARE)

# Select samples
samples <- adni_spare$Sample
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
adni_spare <- left_join(adni_spare,cell_proportions,by="Sample")

# PREPARE FINAL SE
betas_cu <- assays(se)[[1]]
coldata_cu <- as.data.frame(adni_spare)
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

# SAVE
saveRDS(se,'~/Projects/hippoage/data/ADNI/20250314_ADNI_SE.rds')

