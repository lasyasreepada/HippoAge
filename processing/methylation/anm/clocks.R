library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(dnaMethyAge)
library(SummarizedExperiment)

anm <- readRDS('~/Projects/hippoage/data/ANM/anm.rds')

# Extract betas
betas <- as.matrix(assays(anm)[[1]])

# Extract colData 
pheno <- as.data.frame(colData(anm))

# Rename
pheno <- pheno %>%
  dplyr::rename(Age = age.ch1,
         Sample = geo_accession)

# Factor
pheno$Age <- as.integer(pheno$Age)

# Compute clock ages and age accelerations
HorvathS2013 <- methyAge(betas, clock='HorvathS2013', age_info=pheno,simple_mode = TRUE)
HannumG2013 <- methyAge(betas, clock="HannumG2013", age_info=pheno)
HorvathS2018 <- methyAge(betas, clock='HorvathS2018', age_info=pheno,simple_mode = TRUE)
LevineM2018 <- methyAge(betas, clock="LevineM2018", age_info=pheno)
LuA2019 <- methyAge(betas, clock='LuA2019', age_info=pheno)
ShirebyG2020 <- methyAge(betas, clock="ShirebyG2020", age_info=pheno)

library(DunedinPACE)
DunedinPACE <- PACEProjector(betas)

# DunedinPACE
pheno$DunedinPACE <- unlist(DunedinPACE, use.names = FALSE)

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

# Select columns
pheno <- pheno %>%
  select(-c(Sentrix_ID,Sentrix_Position,ID,Deceased,DOD,Education,Race,Ethnicity))

# Write to CSV
write.csv(pheno,'~/Projects/hippoage/data/ANM/clocks.csv')



