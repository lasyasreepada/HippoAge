library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(EnvStats)
library(tableone)
library(stringr)

setwd('/Users/lasyasreepada/Projects/CoMeT/')
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

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
adni_mri <- istaging[istaging$Study=='ADNI',]
adni_mri$SPAREBAG <- adni_mri$SPARE_BA.SVM.RBF. - adni_mri$Age
rm(istaging)
adni_mri <- adni_mri %>%
  select(c(PTID,Date,SPARE_BA.SVM.RBF.,SPAREBAG)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(MRIDATE = Date)

# Format date columns
adni_mri$MRIDATE <- as.Date(adni_mri$MRIDATE)

# Get RID
extract_and_convert <- function(string) {
  last_part <- tail(strsplit(string, "_")[[1]], 1)
  as.integer(last_part)
}

# Example usage
adni_mri$RID <- as.integer(sapply(adni_mri$PTID, extract_and_convert))

# Merge with DNA
adni <- adni %>%
  left_join(adni_mri, by=c('RID')) %>%
  mutate(DNAtoSPARE = abs(interval(SCANDATE, MRIDATE) %>% as.numeric('years'))) %>%
  filter(complete.cases(SPARE_BA.SVM.RBF.)) %>%
  group_by(RID,SCANDATE) %>%
  slice_min(DNAtoSPARE, with_ties = FALSE) %>%
  ungroup()

ad <- adni[adni$DX_nearest_1.0_ordinal==3,]
mci <- adni[adni$DX_nearest_1.0_ordinal==2,]
cu <- adni[adni$DX_nearest_1.0_ordinal==1,]

# SPAREBAG Groups
# Compute mean and SD EAA
mean.sparebag <- mean(adni$SPAREBAG,na.rm=T)
sd.sparebag <- sd(adni$SPAREBAG,na.rm=T)

# Compute confidence intervals
upper.sparebag <- mean.sparebag + 1*sd.sparebag
lower.sparebag <- mean.sparebag - 1*sd.sparebag

# Categorize subjects by epigenetic age group
adni <- adni %>% mutate(GroupSPAREBAG = case_when(SPAREBAG > upper.sparebag ~ "Accelerated",
                                                     SPAREBAG < lower.sparebag ~ "Decelerated",
                                                     (SPAREBAG < upper.sparebag & SPAREBAG > lower.sparebag) ~ "Neutral"))

# Cognitive and SPAREBA correlations
boxplot(sym$PHC_EXF_nearest_CoMeT~sym$GroupSPAREBAG)
boxplot(sym$PHC_MEM_nearest~sym$GroupSPAREBAG)

