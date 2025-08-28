library(dplyr)
library(lme4)
library(RSQLite)
library(data.table)
library(optimx)
library(ggforce)

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

rm(istaging)
adni_mri <- adni_mri %>%
  select(c(PTID,Date,SPARE_BA,Age,Sex,Diagnosis_nearest_2.0)) %>%
  filter(complete.cases(.)) %>%
  dplyr::rename(Age_MRI = Age,
                MRIDATE = Date,
                DX = Diagnosis_nearest_2.0)

adni_mri <- adni_mri[adni_mri$DX=='CN',]

adni_mri <- adni_mri %>%
  group_by(MRIDATE) %>%
  slice_sample(n=1)

# Count the number of timepoints
grouped <- as.data.table(adni_mri)[, count := uniqueN(SPARE_BA), by=PTID]

# Extract longitudinal data with > than 3 timepoints overall
adni_mri <- grouped[count >= 3,]

# LMEM Model
response <- "SPARE_BA"
lmem_formula = as.formula(paste0(response, " ~ Age_MRI + (1 + Age_MRI|PTID)"))

# Fit LMEM
lmem <- lmer(lmem_formula, data=adni_mri, control=lmerControl(optimizer ="Nelder_Mead",
                                                                  calc.derivs = FALSE, optCtrl = list(method = "optimx", starttests = FALSE, kkt = FALSE)))


slopes <- coef(lmem)$PTID['Age_MRI']
names(slopes) <- 'SPARE_BA_slope'
slopes$PTID <- rownames(slopes)
adni_mri <- left_join(adni_mri,slopes,by="PTID")

adni_mri$SPARE_BA_pred <- fitted(lmem)

# Create plot for selected subjects with the individual regression lines overlaid
myplot <- ggplot(adni_mri[1:500,]) + 
  geom_point(aes(x = Age_MRI, y = SPARE_BA), color = "red") +
  geom_smooth(method='lm',aes(x = Age_MRI, y = SPARE_BA_pred)) + 
  facet_wrap_paginate(~PTID, nrow=4, ncol = 4, page=1)

# Calculate the number of pages needed
n_pages <- n_pages(myplot)

# Create a multi-page PDF
pdf("plot.pdf", width = 8.5, height = 11)

for (i in 1:n_pages) {
  print(
    myplot + facet_wrap_paginate(~PTID, nrow=4, ncol = 4, page=i)
  )
}

dev.off()

write.csv(adni_mri,'~/Projects/hippoage/data/ADNI/MRI/20250313_SPAREBAR.csv')

