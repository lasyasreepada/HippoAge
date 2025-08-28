library(dplyr)

# Merge metadata files
mri3meta <- read.csv('~/Projects/hippoage/data/AIBL/MRI/aibl_mri3meta_01-Jun-2018.csv')
mrimeta <- read.csv('~/Projects/hippoage/data/AIBL/MRI/aibl_mrimeta_01-Jun-2018.csv')

# Select necessary columns
mri3meta <- mri3meta %>%
  select(c(RID, SITEID, VISCODE, EXAMDATE))

mrimeta <- mrimeta %>%
  select(c(RID, SITEID, VISCODE, EXAMDATE))

# Concatenate
mrimeta_all <- rbind(mrimeta,mri3meta)
mrimeta_all <- mrimeta_all %>%
  group_by(RID,VISCODE) %>%
  distinct() %>%
  ungroup

write.csv(mrimeta_all,'~/Projects/hippoage/data/AIBL/MRI/MRIMETA_All_4Mar2025.csv',row.names = FALSE)

# Rewrite download column
collections <- read.csv('~/Projects/hippoage/data/AIBL/MRI/AIBL_T1_3_04_2025.csv')
collections <- collections %>%
  select(-c(X.1,X))
collections$Downloaded <- ""

write.csv(collections,'~/Projects/hippoage/data/AIBL/MRI/AIBL_T1_3_04_2025.csv', row.names = FALSE)
