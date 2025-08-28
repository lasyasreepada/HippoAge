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
ashs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/T1_ASHSvols_MTTHK_20250217.csv')
ants <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/thickness_PMTAU_20250218.csv')

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
ashs$MRIDATE <- as.Date(ashs$MRIDATE)
ants$MRIDATE <- as.Date(ants$MRIDATE)

ants <- ants %>%
  select(-c(IMAGEUID.T1,IMAGEUID.T2,IMAGEUID.FLAIR,PHASE_mri,VISCODE_mri,VISCODE2_mri,SEQUENCE_NAME.T1,SEQUENCE_NAME.T2,SEQUENCE_NAME.FLAIR))

# Merge ASHS with ANTS
mri <- merge(ashs, ants, by=c('RID','ID','MRIDATE'))

# Merge with DNA
adni <- adni %>%
  left_join(mri, by=c('RID'), relationship = 'many-to-many') %>%
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

# Select CU only
# adni <- adni[!adni$DIAGNOSIS==1,]

# Remove missing features
adni_temp <- adni %>%
  filter(if_all(c(Left_ITG_inferior_temporal_gyrus_thickness,
                  Left_TMP_temporal_pole_thickness,
                  Left_SPL_superior_parietal_lobule_thickness,
                  Left_PCu_precuneus_thickness,
                  Left_SMG_supramarginal_gyrus_thickness,
                  Left_SFG_superior_frontal_gyrus_thickness,
                  Left_MFG_middle_frontal_gyrus_thickness, 
                  Left_Ent_entorhinal_area_thickness,
                  Left_PHG_parahippocampal_gyrus_thickness,
                  M_BA35_MeanTHK_MSTTHKMT_ASHST1_3T), ~ !is.na(.x)))

# Select cross sectional dataset
adni_temp <- adni_temp %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# Compute CoMeT
adni_temp <- adni_temp %>%
  rowwise() %>%
  mutate(M_Pericalcarine_Thk = mean(c(Left_Calc_calcarine_cortex_thickness,Right_Calc_calcarine_cortex_thickness)),
         M_Insula_Thk = mean(c(Left_PIns_posterior_insula_thickness,Right_PIns_posterior_insula_thickness)),
         M_Cuneus_Thk = mean(c(Left_Cun_cuneus_thickness,Right_Cun_cuneus_thickness)),
         M_Fusiform_Thk = mean(c(Left_FuG_fusiform_gyrus_thickness,Right_FuG_fusiform_gyrus_thickness)),
         M_MedialFrontal_Thk = mean(c(Left_MFC_medial_frontal_cortex_thickness,Right_MFC_medial_frontal_cortex_thickness)),
         M_LateralOccipital_Thk = mean(c(Left_IOG_inferior_occipital_gyrus_thickness,Right_IOG_inferior_occipital_gyrus_thickness,
                                         Left_SOG_superior_occipital_gyrus_thickness,Right_SOG_superior_occipital_gyrus_thickness)),
         M_Precentral_Thk = mean(c(Left_PrG_precentral_gyrus_thickness,Right_PrG_precentral_gyrus_thickness)),
         M_InferiorFrontal_Thk = mean(c(Left_TrIFG_triangular_part_of_the_inferior_frontal_gyrus_thickness,Right_TrIFG_triangular_part_of_the_inferior_frontal_gyrus_thickness,
                                        Left_OpIFG_opercular_part_of_the_inferior_frontal_gyrus_thickness,Right_OpIFG_opercular_part_of_the_inferior_frontal_gyrus_thickness)),
         M_InferiorTemporal_Thk = mean(c(Left_ITG_inferior_temporal_gyrus_thickness,Right_ITG_inferior_temporal_gyrus_thickness)),
         M_TemporalPole_Thk = mean(c(Left_TMP_temporal_pole_thickness,Right_TMP_temporal_pole_thickness)),
         M_SuperiorParietal_Thk = mean(c(Left_SPL_superior_parietal_lobule_thickness,Right_SPL_superior_parietal_lobule_thickness)),
         M_Precuneus_Thk = mean(c(Left_PCu_precuneus_thickness,Right_PCu_precuneus_thickness)),
         M_Supramarginal_Thk = mean(c(Left_SMG_supramarginal_gyrus_thickness,Right_SMG_supramarginal_gyrus_thickness)),
         M_SuperiorFrontal_Thk = mean(c(Left_SFG_superior_frontal_gyrus_thickness,Right_SFG_superior_frontal_gyrus_thickness)),
         M_MidFrontal_Thk = mean(c(Left_MFG_middle_frontal_gyrus_thickness,Right_MFG_middle_frontal_gyrus_thickness)),
         # M_InferiorParietal_Thk = mean(c(ST31TA,ST90TA)),
         M_SuperiorTemporal_Thk = mean(c(Left_STG_superior_temporal_gyrus_thickness,Right_STG_superior_temporal_gyrus_thickness)),
         M_PosteriorCingulate_Thk = mean(c(Left_PCgG_posterior_cingulate_gyrus_thickness,Right_PCgG_posterior_cingulate_gyrus_thickness)),
         M_MiddleTemporal_Thk = mean(c(Left_MTG_middle_temporal_gyrus_thickness,Right_MTG_middle_temporal_gyrus_thickness)),
         M_Parahippocampal_Thk = mean(c(Left_PHG_parahippocampal_gyrus_thickness,Right_PHG_parahippocampal_gyrus_thickness)),
         M_Entorhinal_Thk = sum(c(mean(c(Left_Ent_entorhinal_area_thickness,Right_Ent_entorhinal_area_thickness)),M_BA35_MeanTHK_MSTTHKMT_ASHST1_3T)))
       
# Compute signatures
adni_temp <- adni_temp %>%
  rowwise() %>%
  mutate(CorticalLOAD = mean(c(M_InferiorTemporal_Thk, 
                               M_TemporalPole_Thk, 
                               M_SuperiorParietal_Thk, 
                               M_Precuneus_Thk, 
                               M_Supramarginal_Thk, 
                               M_SuperiorFrontal_Thk,
                               M_MidFrontal_Thk)),
         CorticalEOAD = mean(c(M_Fusiform_Thk,
                               M_SuperiorParietal_Thk,
                               M_Precuneus_Thk,
                               M_SuperiorFrontal_Thk,
                               M_MidFrontal_Thk,
                               M_SuperiorTemporal_Thk,
                               M_PosteriorCingulate_Thk,
                               M_MiddleTemporal_Thk)),
         MTLComposite = mean(c(M_Entorhinal_Thk,
                               M_Parahippocampal_Thk)))

signatures <- c('CorticalLOAD','CorticalEOAD','MTLComposite')

# Define cases and controls based on DX and amyloid status
adni_temp <- adni_temp %>% mutate(Case = case_when(
  (DIAGNOSIS==1) ~ 0, # CU
  (PHC_AMYLOID_STATUS==1 & DIAGNOSIS!=1) ~ 1, # Symptomatic
))

# Loop 
for (sig in signatures) {
  
  # Fit regression to controls
  formula <- as.formula(paste0(sig,'~ Age_MRI + Sex'))
  lmem <- lm(formula, data = subset(adni_temp,Case==0))
  
  # Predict for all subjects
  prediction <- predict(lmem,newdata = adni_temp)
  
  # Adjust values
  residual <- adni_temp[sig] - prediction
  name_residual <- paste0(sig,'_res')
  adni_temp[name_residual] <- as.numeric(unlist(residual))
  
  # Standardize
  zscore <- (residual - mean(unlist(na.omit(adni_temp[adni_temp$Case==0,name_residual])))) / sd(unlist(na.omit(adni_temp[adni_temp$Case==0,name_residual])))
  name_z <- paste0(sig,'_wscore')
  adni_temp[name_z] <- as.numeric(unlist(zscore))
  
}

# COMPUTE CoMeT
signatures.w <- grep("OAD_wscore$", names(adni_temp), value=TRUE)
for (sig in signatures.w) {
  comet <- adni_temp[sig] - adni_temp['MTLComposite_wscore']
  name_comet <- paste0(gsub('_wscore','',sig),'_CoMeT')
  adni_temp[name_comet] <- comet
}

# Final Data
adni_comet <- adni_temp[!is.na(adni_temp$Case),]

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
saveRDS(se,'~/Projects/hippoage/data/ADNI/20250316_COMET_SE.rds')
