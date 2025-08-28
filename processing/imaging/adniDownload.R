library(dplyr)
library(data.table)
library(SummarizedExperiment)
library(lubridate)

# READ DATA
adni_dna <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/Methylation Data - Anil + Lasya/adni.rds')
adni_mri <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/T1_ASHSvols_MTTHK_20241108.csv')
dx <- read.csv('~/Projects/hippoage/data/ADNI/DXSUM_21Jan2025.csv')

# GET COLDATA
adni_dnam <- as.data.frame(colData(adni_dna))
rownames(adni_dnam) <- NULL

# TRANSFORM COLUMNS
adni_dnam$DateDrawn <- as.Date(adni_dnam$DateDrawn)
adni_mri$MRIDATE <- as.Date(adni_mri$MRIDATE)
dx$EXAMDATE <- as.Date(dx$EXAMDATE)

# MERGE DX
dx <- dx %>%
  select(c(RID,PTID,VISCODE2,EXAMDATE,DIAGNOSIS))

adni <- adni_dnam %>%
  left_join(dx, by = c("RID"),relationship = 'many-to-many') %>%
  mutate(DNAtoDX = abs(interval(DateDrawn, EXAMDATE) %>% as.numeric('years'))) %>% # Calculate the absolute date difference
  group_by(RID, DateDrawn) %>% # Group by ID and Date1
  slice_min(DNAtoDX, with_ties = FALSE) %>%  # Select the row with the smallest date difference
  ungroup()

# MERGE MRI
adni_temp <- adni %>%
  left_join(adni_mri, by=c('RID'), relationship='many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DateDrawn, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# RECOVER MISSING MRI  
missing <- unique(adni_temp[(is.na(adni_temp$ID) & (adni_temp$DIAGNOSIS==1)),c('PTID','RID','DateDrawn','EXAMDATE','DIAGNOSIS','DNAtoMRI')])
sus <- adni_temp[(!is.na(adni_temp$DNAtoMRI) & adni_temp$DNAtoMRI > 3 & (adni_temp$DIAGNOSIS==1)),c('PTID','RID','DateDrawn','EXAMDATE','DIAGNOSIS','DNAtoMRI')]

# Missing and mismatched individuals
to_download <- rbind(missing,sus)
write.csv(to_download,'~/Projects/hippoage/data/ADNI/to_download.csv')

# Get a list of unique IDs
ids <- as.list(unique(to_download$PTID))
fwrite(ids, "~/Projects/hippoage/data/ADNI/to_download_ids.csv", col.names = FALSE)

# IDA SEARCH
search <- read.csv('~/Downloads/idaSearch_2_15_2025-2.csv')

# Include only MPRAGE and SPGR
search_filtered <- search[grepl("RAGE|SPGR", search$Description),]

# Select higher UID for repeat images
search_filtered <- search_filtered %>%
  group_by(Subject.ID,Study.Date) %>%
  slice_max(Image.ID, with_ties = FALSE) %>%
  ungroup()

# REFORMAT
library(tidyr)

# Split Image Protocol
search_filtered <- search_filtered %>%
  separate(Imaging.Protocol, into = c("Field_Strength", "Weighting"), sep = ";") %>%
  mutate(Field_Strength = sub("Field Strength=", "", Field_Strength),
         Weighting = sub("Weighting=", "", Weighting))

# Get RID
search_filtered$RID <- as.integer(sub(".*_", "", search_filtered$Subject.ID))

# Convert to Date
search_filtered$MRIDATE <- as.Date(search_filtered$Study.Date, format="%m/%d/%Y")

# Select necessary columns
search_filtered <- search_filtered %>%
  select(c(RID,MRIDATE,Image.ID))

# MERGE MRI SEARCH
adni_temp_long <- to_download %>%
  left_join(search_filtered, by=c('RID'), relationship='many-to-many') %>%
  mutate(DNAtoMRI = abs(interval(DateDrawn, MRIDATE) %>% as.numeric('years'))) %>%
  group_by(RID, DateDrawn) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

adni_temp_x <- adni_temp2 %>%
  group_by(RID) %>%
  slice_min(DNAtoMRI, with_ties = FALSE) %>%
  ungroup()

# Get list of UIDs to download
uids <- as.list(unique(adni_temp_x$Image.ID))
fwrite(uids, "~/Projects/hippoage/data/ADNI/to_download_uids.csv", col.names = FALSE)

# MERGE METADATA
library(dplyr)

# Merge metadata files
mri3meta <- read.csv('~/Projects/hippoage/data/ADNI/MRI/MR_Image_Acquisition/MRI3META_16Feb2025.csv')
mrimeta <- read.csv('~/Projects/hippoage/data/ADNI/MRI/MR_Image_Acquisition/MRIMETA_16Feb2025.csv')

# Select necessary columns
mri3meta <- mri3meta %>%
  select(c(PHASE, PTID, RID, VISCODE, VISCODE2, FIELD_STRENGTH, EXAMDATE))

mrimeta <- mrimeta %>%
  select(c(PHASE, PTID, RID, VISCODE, VISCODE2, FIELD_STRENGTH, EXAMDATE))

# Concatenate
mrimeta_all <- rbind(mrimeta,mri3meta)

write.csv(mrimeta_all,'~/Projects/hippoage/data/ADNI/MRI/MR_Image_Acquisition/MRIMETA_All_16Feb2025.csv')


