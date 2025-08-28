library(dplyr)
library(tableone)

# READ DATA
abc <- read.csv('~/Projects/hippoage/data/INDD/20250530_abc_pheno.csv')
adni <- read.csv('~/Projects/hippoage/data/ADNI/20250427_adni_pheno.csv')

# ABC
abc_vars <- c('Age')
abc_cat <- c()

# Remove empty strings
abc$APOE <- ifelse(abc$APOE == "", NA, abc$APOE)

# Create table 1
abc_t1_dx <- CreateTableOne(vars = abc_vars, strata = 'Dx', data = abc, factorVars = abc_cat)
abc_vars <- c(abc_vars,'Dx')
abc_cat <- c(abc_cat,'Dx')
abc_t1 <- CreateTableOne(vars = abc_vars, data = abc, factorVars = abc_cat)

# ADNI
adni_vars <- c('Age_DNA','DNAtoMRI','')
adni_cat <- c('Sex','APOE_Genotype','APOE4_Alleles')

adni$DX_nearest_1.0 <- factor(adni$DX_nearest_1.0, levels = c('CN','MCI','Dementia')) 

# Create table 1
adni_t1_dx <- CreateTableOne(vars = adni_vars, strata = 'DX_nearest_1.0', data = adni, factorVars = adni_cat)
adni_t1 <- CreateTableOne(vars = adni_vars, data = adni, factorVars = adni_cat)


# AIBL
# Characteristics
aibl_vars <- c('Age','Sex','apoe','Demographic.Years.of.Education','Neuropsych.CDR.SOB')
aibl_cat <- c('Sex','apoe','Demographic.Years.of.Education')

# FACTOR
aibl$Diagnosis <- factor(aibl$Diagnosis, levels = c('HC','MCI','AD')) 
aibl$Demographic.Years.of.Education <- factor(aibl$Demographic.Years.of.Education, levels = c('0-6','7-8','9-12','13-15','15+')) 

# Create table 1
aibl_t1_dx <- CreateTableOne(vars = aibl_vars, strata = 'Diagnosis', data = aibl, factorVars = aibl_cat)
aibl_t1 <- CreateTableOne(vars = vars, data = abc, factorVars = cats)


