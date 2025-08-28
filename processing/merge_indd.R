library(dplyr)
library(readxl)
library(ggplot2)
library(tableone)
library(ggVennDiagram)
library(lubridate)

# READ DATA
se <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/indd.rds')
dna <- read_excel('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/MethylationMasterIDs_14Jan2025.xlsx')
mri <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC_MRI3T_ASHST1T2_measurements_20241220.csv')

# GET COLDATA
indd <- as.data.frame(colData(se))
rownames(indd) <- NULL

indd <- indd %>% 
  select(-c(Sentrix_ID,Sentrix_Position)) %>%
  dplyr::rename(Methylome1_DataID = Sample)

# Separate Missing
indd_nona <- indd[!is.na(indd$ID),]
indd_na <- indd[is.na(indd$ID),]
  
# Pull data
indd_na <- left_join(indd_na,dna,by='Methylome1_DataID')
indd_na$ID <- indd_na$Methylome1_SampleID
indd_na$DNAdate <- as.Date(indd_na$Methylome1_SampleDate, format = "%Y/%m/%d")
indd_na <- indd_na %>%
  select(Methylome1_DataID,ID,Sample_Well,Sample_Plate,DNAdate,Version)

indd_nona$DNAdate <- as.Date(indd_nona$DNAdate, format = "%m/%d/%Y")

# CONCAT
indd <- rbind(indd_nona,indd_na)
indd$INDDID <- (sapply(strsplit(indd$ID, "_"), `[`, 1))

# MERGE DEMOGRAPHICS
dem <- read_excel('/Users/lasyasreepada/Projects/hippoage/data/INDD/INDD_Methylation_Demographics_DX.xlsx')
dem$INDDID <- as.character(dem$INDDID)
indd <- left_join(indd, dem, by='INDDID')

# Compute Age
indd$Age <- time_length(interval(as.Date(indd$DOB),as.Date(indd$DNAdate),), unit = "years")
  
# SAVE
write.csv(indd,'/Users/lasyasreepada/Projects/hippoage/data/INDD/temp.csv')
  
# MRI
mri_ids <- unique(mri$INDD)
mri_ids <- unique(gsub("_", ".", mri_ids))

# SETS

# Intersection
length(intersect(indd_ids,mri_ids))

# Difference
list <- setdiff(ids_indd_base,intersect(ids_DNA,ids_indd_base))
write.csv(list,'~/Projects/hippoage/data/indd_DNAm_no_MRI.csv')

# VENN
venn_data <- list(New = ids_DNA_base, Old = ids_indd_base)

# Plot the Venn diagram
ggVennDiagram(venn_data) +
  scale_fill_gradient(
    low = "lightblue",  # Color for low values
    high = "darkblue",  # Color for high values
  ) +
  ggtitle('Venn Diagram of Imaging and Methylation overlap')

# CHARACTERISTICS
myVars <- c('AgeatMRI','Sex','Race','Ethnicity','APOE','Education','CDRSum','M_Hippo_MeanTHK_MSTTHKMT_ASHST1_3T',"M_ERC_MeanTHK_MSTTHKMT_ASHST1_3T","M_PHC_MeanTHK_MSTTHKMT_ASHST1_3T")
catVars <- c('Sex','Race','Ethnicity','APOE')

tab1_dx <- CreateTableOne(vars = myVars, strata = 'Dx', data = indd_3T, factorVars = catVars)


