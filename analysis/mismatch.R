library(dplyr)

adni <- readRDS('~/Projects/hippoage/data/ADNI/ADNI_SPAREBA_SE.rds')

coldata <- as.data.frame(colData(adni))
rownames(coldata) <- NULL

coldata$SPAREBAG <- coldata$SPARE_BA.SVM.RBF. - coldata$Age_MRI

# # Categorize subjects by biological age groups
# coldata <- coldata %>% mutate(GroupLevine2018 = case_when(LevineM2018_res > 0 ~ "Accelerated",
#                                                      LevineM2018_res < 0 ~ "Decelerated",
#                                                      (LevineM2018_res == 0) ~ "Neutral"))
# Categorize subjects by biological age groups
# DUNEDIN
# Compute mean and SD EAA
mean.dunedin <- mean(coldata$DunedinPACE,na.rm=T)
sd.dunedin <- sd(coldata$DunedinPACE,na.rm=T)

# Compute confidence intervals
upper.dunedin <- mean.dunedin + 1*sd.dunedin
lower.dunedin <- mean.dunedin - 1*sd.dunedin

# Categorize subjects by epigenetic age group
coldata <- coldata %>% mutate(GroupDunedinPACE = case_when(DunedinPACE > upper.dunedin ~ "Accelerated",
                                                     DunedinPACE < lower.dunedin ~ "Decelerated",
                                                     (DunedinPACE < upper.dunedin & DunedinPACE > lower.dunedin) ~ "Neutral"))

# SPAREBA
# Compute mean and SD EAA
mean.spareba <- mean(coldata$SPAREBAG,na.rm=T)
sd.spareba <- sd(coldata$SPAREBAG,na.rm=T)

# Compute confidence intervals
upper.spareba <- mean.spareba + 1*sd.spareba
lower.spareba <- mean.spareba - 1*sd.spareba

# Categorize subjects by brain age groups
coldata <- coldata %>% mutate(GroupSPAREBA = case_when(SPAREBAG > upper.spareba ~ "Accelerated",
                                                       SPAREBAG < lower.spareba ~ "Decelerated",
                                                       (SPAREBAG < upper.spareba & SPAREBAG > lower.spareba) ~ "Neutral"))
# Table
cu <- coldata[coldata$DIAGNOSIS==1,]
table(cu$GroupDunedinPACE,cu$GroupSPAREBA)
cor(cu$DunedinPACE,cu$SPAREBAG)

