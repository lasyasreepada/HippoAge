library(dplyr)
library(readxl)
library(SummarizedExperiment)

# READ 
rush.long <- ''
radc.rosmap <- read.csv("/Users/lasyasreepada/Library/CloudStorage/Box-Box/Rush_data/Rush_Jan2025_shared_files/RADC_MRI_ROS-MAP_share_2024-12-03/RADC_MRI_ROS-MAP_share_2024-12-03.csv")
rosmap <- readRDS("/Users/lasyasreepada/Library/CloudStorage/Box-Box/Data/20240613_ROSMAP_SE_New.rds")
radc.long <- read_excel("/Users/lasyasreepada/Library/CloudStorage/Box-Box/Rush_data/Rush_Jan2025_shared_files/dataset_1272_long_01-09-2025.xlsx")
  
# GET PHENO
pheno <- as.data.frame(colData(rosmap))
rownames(pheno) <- NULL

# GET IDS
ids.mri <- unique(as.integer(sub("\\'", "", radc.rosmap$projid)))
ids.radc <- unique(as.integer(radc.long$projid))
ids.dna <- unique(pheno$subject_id)

# VENN
length(intersect(ids.dna,ids.mri))
length(intersect(ids.dna,ids.mri.radc))