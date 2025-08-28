library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(dnaMethyAge)
library(SummarizedExperiment)

indd <- readRDS('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/indd.rds')

# Remove rows and features with too much missingness
cleanMatrixForClusterSE <- function(se, f_row = 0.5, f_col = 0.5) {
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

indd <- cleanMatrixForClusterSE(indd,f_row = 1)
pheno <- as.data.frame(colData(indd))

# Read new annotation file
annotation <- read_excel('/Users/lasyasreepada/Library/CloudStorage/Box-Box/ABC/Methylation/MethylationMasterIDs_14Jan2025.xlsx')

# Identify missing
missing <- pheno[is.na(pheno$ID),'Sample']
recover <- missing %in% annotation$Methylome1_DataID

# Standardize INDDID
pheno$INDDID <- (sapply(strsplit(pheno$ID, "_"), `[`, 1))
annotation$INDDID <- as.character(annotation$INDDID)

# Standardize Date
pheno$DNAdate <- as.Date(pheno$DNAdate,format = "%m/%d/%Y")
annotation$Methylome1_SampleDate <- as.Date(annotation$Methylome1_SampleDate)

lookup <- annotation[annotation$Methylome1_DataID %in% missing,]
lookup <- lookup %>%
  dplyr::rename(Sample = Methylome1_DataID,
                DNAdate = Methylome1_SampleDate,
                ID = Methylome1_SampleID) %>%
  select(c(Sample,ID,INDDID,DNAdate))

common_cols <- intersect(names(pheno), names(lookup))

# Fill NAs dynamically for only common columns
pheno_filled <- pheno %>%
  left_join(lookup, by = "Sample", suffix = c("", "_lookup")) %>%  # Join lookup table
  mutate(across(all_of(common_cols), ~ ifelse(is.na(.), get(paste0(cur_column(), "_lookup")), .))) %>% # Fill NAs
  select(-ends_with("_lookup"))  # Remove temporary lookup columns

pheno_filled$DNAdate <- as.Date(pheno_filled$DNAdate)  
rownames(pheno_filled) <- pheno_filled$Sample

# ADD DEMOGRAPHICS
dem <- read_excel('~/Projects/hippoage/data/INDD/Demographics.xlsx')
dem$INDDID <- as.character(dem$INDDID)
pheno_filled <- left_join(pheno_filled,dem,by='INDDID')

# COMPUTE AGE
library(lubridate)
pheno_filled$Age <- interval(pheno_filled$DOB, pheno_filled$DNAdate) %/% years(1)

# CREATE NEW SE
betas <- assays(indd)[[1]]
coldata <- pheno_filled
qc <- metadata(indd)
se <- SummarizedExperiment(assays = betas, colData=coldata, metadata = qc)

# SAVE  
saveRDS(se, "/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/indd.rds")

# Compute Clocks
indd <- readRDS("/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/indd.rds")
library(dnaMethyAge)

# Extract betas
betas <- assays(indd)[[1]]

# Extract colData 
pheno <- as.data.frame(colData(indd))

# Select individuals not missing age
pheno <- pheno[!is.na(pheno$Age),]
betas <- betas[,rownames(pheno)]

# Compute clock ages and age accelerations
HorvathS2013 <- methyAge(betas, clock='HorvathS2013', age_info=pheno,simple_mode = TRUE)
HannumG2013 <- methyAge(betas, clock="HannumG2013", age_info=pheno)
HorvathS2018 <- methyAge(betas, clock='HorvathS2018', age_info=pheno,simple_mode = TRUE)
LevineM2018 <- methyAge(betas, clock="LevineM2018", age_info=pheno)
LuA2019 <- methyAge(betas, clock='LuA2019', age_info=pheno)
ShirebyG2020 <- methyAge(betas, clock="ShirebyG2020", age_info=pheno)

library(DunedinPACE)
DunedinPACE <- PACEProjector(betas)

plot(pheno$Age,pheno$DunedinPACE,
     xlab = "Age", ylab = "DunedinPACE", 
     main = "DunedinPACE")
abline(h = 1, col = "red", lty = 2)

# CONCATENATE
# Horvath 1
pheno$HorvathS2013 <- HorvathS2013$mAge
pheno$HorvathS2013_res <- HorvathS2013$Age_Acceleration

# Hannum
pheno$HannumG2013 <- HannumG2013$mAge
pheno$HannumG2013_res <- HannumG2013$Age_Acceleration

# Horvath 2
pheno$HorvathS2018 <- HorvathS2018$mAge
pheno$HorvathS2018_res <- HorvathS2018$Age_Acceleration

# Levine (Pheno)
pheno$LevineM2018 <- LevineM2018$mAge
pheno$LevineM2018_res <- LevineM2018$Age_Acceleration

# Lu (Telomere)
pheno$LuA2019 <- LuA2019$mAge
pheno$LuA2019_res <- LuA2019$Age_Acceleration

# Shireby (Cortical)
pheno$ShirebyG2020 <- ShirebyG2020$mAge
pheno$ShirebyG2020_res <- ShirebyG2020$Age_Acceleration

# DunedinPACE
pheno$DunedinPACE <- unlist(DunedinPACE, use.names = FALSE)

# Select columns
pheno <- pheno %>%
  select(-c(Sentrix_ID,Sentrix_Position,ID,Deceased,DOD,Education,Race,Ethnicity))

# Write to CSV
write.csv(pheno,'~/Projects/hippoage/data/INDD/indd_clocks_20250201.csv')



