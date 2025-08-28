library(dplyr)
library(tableone)

# ADNI 
adni <- readRDS('~/Projects/hippoage/data/ADNI/20250312_ADNI_SE.rds')
coldata_adni <- as.data.frame(colData(adni))
coldata_adni$PHC_AMYLOID_STATUS <- as.factor(coldata_adni$PHC_AMYLOID_STATUS)
coldata_adni$APOE4_Alleles <- as.factor(coldata_adni$APOE4_Alleles)

vars_adni <- c('Age_DNA','DNAtoMRI','Sex','APOE_Genotype','APOE4_Alleles','Education_Years','Race','Ethnicity','CDRSB','PHC_AMYLOID_STATUS','M_Hippo_VOL_ASHST1.combat')
cat_adni <- c('Sex','Race','Ethnicity','APOE_Genotype','APOE4_Alleles','PHC_AMYLOID_STATUS')
medians_adni <- c('CDRSB','Education_Years')
table1_adni <- CreateTableOne(vars = vars_adni, data = coldata_adni, factorVars = cat_adni)

# ABC 
abc <- readRDS('~/Projects/hippoage/data/INDD/20250312_ABC_SE.rds')
coldata_abc <- as.data.frame(colData(abc))
coldata_abc$APOE4_Alleles <- as.factor(coldata_abc$APOE4_Alleles)
vars_abc <- c('Age_DNA','DNAtoMRI','Sex','Race','Ethnicity','APOE_Genotype','APOE4_Alleles','Education','CDRSum..StandardCDRSum.','Amy_Status','M_Hippo_VOL_ASHST1')
cat_abc <- c('Sex','Race','Ethnicity','APOE_Genotype','APOE4_Alleles','Amy_Status')
medians_abc <- c('CDRSum..StandardCDRSum.','Education')
table1_abc <- CreateTableOne(vars = vars_abc, data = coldata_abc, factorVars = cat_abc)

# AIBL
aibl <- read.csv('~/Projects/hippoage/data/AIBL/20250315_aibl_pheno.csv')
vars_aibl <- c('Age_DNA','DNAtoMRI','Sex','Demographic.ApoE.genotype','Demographic.Years.of.Education','Neuropsych.CDR.SOB','M_Hippo_VOL_ASHST1')
cat_aibl <- c('Sex','Demographic.ApoE.genotype','Demographic.Years.of.Education')
medians_aibl <- c('Neuropsych.CDR.SOB')
table1_aibl <- CreateTableOne(vars = vars_aibl, data = aibl, factorVars = cat_aibl)

# PRINT
print(table1_adni, nonnormal = medians_adni)
print(table1_abc, nonnormal = medians_abc)
print(table1_aibl, nonnormal = medians_aibl)



