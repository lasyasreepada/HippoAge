library(dplyr)

# READ DATA
fs <- readRDS('~/Projects/hippoage/data/ADNI/20250319_COMET_SE.rds')
ants <- readRDS('~/Projects/hippoage/data/ADNI/20250316_COMET_SE.rds')

# Get coldata
# coldata_fs <- as.data.frame(colData(fs))
coldata_fs <- as.data.frame(colData(fs))
coldata_ants <- as.data.frame(colData(ants))

# Select CoMeT columns
comet_fs <- coldata_fs %>%
  select(c(RID,Case,CorticalLOAD_CoMeT,CorticalEOAD_CoMeT)) %>%
  dplyr::rename(CorticalLOAD_CoMeT_FS = CorticalLOAD_CoMeT,
                CorticalEOAD_CoMeT_FS = CorticalEOAD_CoMeT)

comet_ants <- coldata_ants %>%
  select(c(RID,Case,CorticalLOAD_CoMeT,CorticalEOAD_CoMeT)) %>%
  dplyr::rename(CorticalLOAD_CoMeT_ANTS = CorticalLOAD_CoMeT,
                CorticalEOAD_CoMeT_ANTS = CorticalEOAD_CoMeT)

# Merge
comet <- merge(comet_fs,comet_ants,by='RID',all=TRUE)

# Select rows without NA
comet <- comet %>%
  filter(if_all(c(CorticalLOAD_CoMeT_FS,CorticalEOAD_CoMeT_FS,CorticalLOAD_CoMeT_ANTS,CorticalEOAD_CoMeT_ANTS), ~ !is.na(.x)))

# Plots
# CorticalLOAD_CoMeT
cor.test(comet$CorticalLOAD_CoMeT_FS,comet$CorticalLOAD_CoMeT_ANTS)
plot(comet$CorticalLOAD_CoMeT_FS,comet$CorticalLOAD_CoMeT_ANTS)

# CorticalEOAD_CoMeT
cor.test(comet$CorticalEOAD_CoMeT_FS,comet$CorticalEOAD_CoMeT_ANTS)
plot(comet$CorticalEOAD_CoMeT_FS,comet$CorticalEOAD_CoMeT_ANTS)

cases_fs <- coldata_fs[coldata_fs$Case==1,]
cases_ants <- coldata_ants[coldata_ants$Case==0,]
