library(RSQLite)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(invgamma)
library(lmtest)
library(readxl)

# Setup SQLite driver
sqlite.driver <- dbDriver("SQLite")

# Establish database connection
file <- "~/Projects/predict4ad/data/istaging.db"
db <- dbConnect(sqlite.driver, dbname=file)

# Read in brain age data
dbListTables(db)
data <- dbReadTable(db,"istaging")

# Disconnect SQLite
dbDisconnect(db)

# Read in epigenetic clock data
adni <- read.csv('~/Projects/havana/data/ADNI/ADNI_DNAM_MRI_COG_LONG_H4.csv')
# aibl <- read_excel('/Users/lasyasreepada/Library/CloudStorage/Box-Box/AIBL/Corey_age_project_DNAmAge_data_cognition_imaging_20240916.xlsx')

# Complete cases
data.1 <- data[complete.cases(data[,c('PTID','Age','Visit_Code','SPARE_BA','SPARE_BA.SVM.RBF.','Study')]),]

# Compute SPARE BAG
data.1 <- data.1 %>% rowwise() %>% mutate(SPARE_BAG = SPARE_BA.SVM.RBF.- Age)

# AIBL
# Clean
aibl <- subset(aibl, select = -c(...31,...32,AIBL_ID...33)) 
aibl <- aibl %>% dplyr::rename(PTID := AIBL_ID...1)

data.2 <- data.1[data.1$Study=='AIBL',c('PTID','Age','Visit_Code','SPARE_BA.SVM.RBF.','SPARE_BAG')]
data.2$PTID <- as.double(data.2$PTID)

# Merge
data.aibl <- left_join(subset(aibl,select=-Age), data.2, by=c('PTID'))
data.aibl.x <- data.aibl %>%
  group_by(PTID) %>%
  dplyr::slice(which.min(Age)) %>% # For subjects with DNAm, select row with min Years_to_DNAm
  ungroup()

# ADNI
# Complete Cases
adni.1 <- adni[complete.cases(adni[,c('RID','VISCODE2','AGE_SAMPLE_nearest','DX_nearest_1.0_ordinal')]),]
adni.1$DX_nearest_1.0 = factor(adni.1$DX_nearest_1.0,levels=c('CN','MCI','Dementia'))

# Derive RID from PTID
data.2 <- data.1[data.1$Study=='ADNI',c('PTID','Age','Visit_Code','SPARE_BA','SPARE_BA.SVM.RBF.','SPARE_BAG')]
data.2$RID <- as.integer(sapply(data.2$PTID, function(x) strsplit(x, "_")[[1]][3], USE.NAMES=FALSE))

# Merge brain age data with epigenetic data
data.2 <- data.2 %>% dplyr::rename(VISCODE2 = Visit_Code)
data.3 <- left_join(adni.1, data.2, by=c('RID','VISCODE2'),relationship = 'many-to-many')
data.adni <- data.3[complete.cases(data.3[,c('RID','SPARE_BA.SVM.RBF.','AgeHorvath2018_nearest')]),]

data.adni.x <- data.adni %>%
  group_by(RID) %>%
  dplyr::slice(which.min(Years_to_DNAm)) %>% # For subjects with DNAm, select row with min Years_to_DNAm
  ungroup()
  
# Pool
data.adni.x.temp <- subset(data.adni.x,select=c("RID","Age","PTGENDER_ordinal","DX_nearest_1.0_ordinal","AMYLOID_SUMMARY","AgeHorvath2013_nearest","AgeShireby2020_nearest","AgeHorvath2018_nearest", "DunedinPACE_nearest","AgeHannum2013_nearest","AgeAccHorvath2013_nearest","AgeAccShireby2020_nearest","AgeAccHorvath2018_nearest","AgeAccHannum2013_nearest","SPARE_BA.SVM.RBF.",'SPARE_BAG',"CDRSB"))
data.aibl.x.temp <- subset(data.aibl.x,select=c('PTID','Age','Sex','Diagnosis','Centiloid','HorvathS2013','ShirebyG2020','HorvathS2018','DunedinPACE','HannumG2013','HorvathS2013_res','ShirebyG2020_res','HorvathS2018_res','HannumG2013_res','SPARE_BA.SVM.RBF.','SPARE_BAG','Neuropsych.CDR.SOB'))

data.adni.x.temp$RID <- paste(data.adni.x.temp$RID, "ADNI", sep="_")
data.aibl.x.temp$PTID <- paste(data.aibl.x.temp$PTID, "AIBL", sep="_")

data.adni.x.temp$Study <- "ADNI"
data.aibl.x.temp$Study <- "AIBL"

colnames <- c('PTID','Age','Sex','DX','Amyloid','HorvathS2013','ShirebyG2020','HorvathS2018','DunedinPACE','HannumG2013','HorvathS2013_res','ShirebyG2020_res','HorvathS2018_res','HannumG2013_res','SPAREBA2','SPAREBAG','CDRSOB','Study')

names(data.adni.x.temp) <- colnames
names(data.aibl.x.temp) <- colnames

data.pooled <- rbind(data.adni.x.temp,data.aibl.x.temp)

# Pooled Plots

# Correlation
ggplot(data.pooled,aes(ShirebyG2020,SPAREBA2,color=Study)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BA vs Shireby", x='Shireby', y='SPARE BA') + 
  geom_vline(xintercept = mean(data.pooled$ShirebyG2020)) + 
  geom_hline(yintercept = mean(data.pooled$SPAREBA2))

# Bland Altman
library(blandr)
blandr.statistics(data.pooled$HorvathS2018, data.pooled$SPAREBA2, sig.level=0.95)



# Plots
ggplot(data.aibl,aes(Age,SPARE_BA.SVM.RBF.,group = PTID,color = Diagnosis)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="SPARE BA vs Age - Longitudinal", x='Age',y='SPARE BA')

ggplot(data.aibl.x,aes(Age,SPARE_BA.SVM.RBF.,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BA vs Age - Cross-Sectional", x='Age',y='SPARE BA')

ggplot(data.aibl,aes(Age,SPARE_BAG,group=PTID,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="SPARE BAG vs Age - Longitudinal",x='Age',y='SPARE BAG')

ggplot(data.aibl.x,aes(Age,SPARE_BAG,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BAG vs Age - Cross-Sectional",x='Age',y='SPARE BAG')

ggplot(data.aibl.x,aes(Age,HorvathS2018,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + labs(title="Horvath 2 vs Age",x='Age',y='Horvath 2')

ggplot(data.aibl.x,aes(Age,HorvathS2018_res,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + labs(title="Horvath 2 Age Gap",x='Age',y='Horvath 2 Age Gap')

ggplot(data.aibl.x,aes(HorvathS2018,SPARE_BA.SVM.RBF.,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BA vs Horvath 2", x='Horvath 2', y='SPARE BA') + 
  geom_vline(xintercept = mean(data.aibl.x$HorvathS2018)) + 
  geom_hline(yintercept = mean(data.aibl.x$SPARE_BA.SVM.RBF.))

ggplot(data.aibl.x,aes(HorvathS2018_res,SPARE_BAG,colour = Diagnosis)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BAG vs Horvath 2 Age Gap",x='Horvath 2 Age Gap',y='SPARE BAG') + 
  geom_vline(xintercept = mean(data.aibl.x$HorvathS2018_res)) + 
  geom_hline(yintercept = mean(data.aibl.x$SPARE_BAG))

ggplot(data.adni,aes(Age,SPARE_BA.SVM.RBF.,group = RID,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="SPARE BA vs Age - Longitudinal",x='Age',y='SPARE BA')

ggplot(data.adni.x,aes(Age,SPARE_BA.SVM.RBF.,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BA vs Age - Cross-Sectional", x='Age',y='SPARE BA')

ggplot(data.adni,aes(Age,SPARE_BAG,group=RID,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="SPARE BAG vs Age - Longitudinal",x='Age',y='SPARE BAG')

ggplot(data.adni.x,aes(Age,SPARE_BAG,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BAG vs Age - Cross-Sectional",x='Age',y='SPARE BAG')

ggplot(data.adni,aes(AGE_SAMPLE_nearest,AgeHorvath2018_nearest,group = RID,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="Horvath 2 - Longitudinal",x='Age',y='Horvath 2')

ggplot(data.adni.x,aes(AGE_SAMPLE_nearest,AgeHorvath2018_nearest,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + labs(title="Horvath 2 vs Age - Cross Sectional",x='Age',y='Horvath 2')

ggplot(data.adni,aes(AGE_SAMPLE_nearest,AgeAccHorvath2018_nearest,group = RID,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="Horvath 2 Age Gap - Longitudinal",x='Age',y='Horvath 2')

ggplot(data.adni.x,aes(AGE_SAMPLE_nearest,AgeAccHorvath2018_nearest,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + labs(title="Horvath 2 Age Gap - Cross Sectional",x='Age',y='Horvath 2')

ggplot(data.adni,aes(AgeHorvath2018_nearest,SPARE_BA.SVM.RBF.,group = RID,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE) + labs(title="SPARE BA vs Horvath 2 - Longitudinal", x='Horvath 2', y='SPARE BA')

ggplot(data.adni.x,aes(AgeHorvath2018_nearest,SPARE_BA.SVM.RBF.,group = RID,colour = DX_nearest_1.0)) +
  geom_point(show.legend = FALSE) + labs(title="SPARE BA vs Horvath 2 - Cross-Sectional",x='Horvath 2',y='SPARE BA') + geom_vline(xintercept = mean(data.adni.x$AgeHorvath2018_nearest)) + geom_hline(yintercept = mean(data.adni.x$SPARE_BA.SVM.RBF.))

ggplot(data.adni.x,aes(DX_nearest_1.0,SPARE_BA.SVM.RBF.,colour = DX_nearest_1.0)) +
  geom_boxplot(show.legend = FALSE) + labs(title="SPARE BA vs Clinical DX", x='Clinical DX', y='SPARE BA')

ggplot(data.adni.x,aes(DX_nearest_1.0,AgeHorvath2018_nearest,colour = DX_nearest_1.0)) +
  geom_boxplot(show.legend = FALSE) + labs(title="Horvath 2 vs Clinical DX", x='Clinical DX', y='Horvath 2')

ggplot(data.adni.x,aes(DX_nearest_1.0,SPARE_BAG,colour = DX_nearest_1.0)) +
  geom_boxplot(show.legend = FALSE) + labs(title="SPARE BAG vs Clinical DX", x='Clinical DX', y='SPARE BAG')

ggplot(data.adni.x,aes(DX_nearest_1.0,AgeAccHorvath2018_nearest,colour = DX_nearest_1.0)) +
  geom_boxplot(show.legend = FALSE) + labs(title="Horvath 2 Age Gap vs Clinical DX", x='Clinical DX', y='Horvath 2 Age Gap')

# Statistical Analyses
anova <- aov(SPARE_BA.SVM.RBF. ~ DX_nearest_1.0 + Age + PTGENDER, data=adni.5)

y <- adni.5$SPARE_BA.SVM.RBF.
x <- adni.5$Age

# 1. Fit OLS model and plot residuals
model <- lm(y ~ x, data = adni.5)
plot(model$fitted.values, residuals(model), xlab = "Fitted values", ylab = "Residuals")
abline(h = mean(y), col = "red")
summary(model)$r.squared # R2 = 0.0009691044

# 2. Breusch-Pagan test
bptest(model)

# 3. White's test
bptest(model, ~ fitted.values(model) + I(fitted.values(model)^2))

# 4. Goldfeld-Quandt test
gqtest(model, order.by = ~ x, data = adni.5)

# # 5. Levene's test
# install.packages("car")
# library(car)
# leveneTest(y ~ factor(x_group), data = your_data)

# # 6. Weighted Least Squares (WLS)
# weights <- 1 / fitted(model)^2
# wls_model <- lm(y ~ x, data = your_data, weights = weights)
# summary(wls_model)

y <- adni.5$SPARE_BAG
x <- adni.5$Age

# 1. Fit OLS model and plot residuals
model <- lm(y ~ x, data = adni.5)
plot(model$fitted.values, residuals(model), xlab = "Fitted values", ylab = "Residuals")
abline(h = mean(y), col = "red")
summary(model)$r.squared # R2 = 0.003840183

# 2. Breusch-Pagan test
bptest(model)

# 3. White's test
bptest(model, ~ fitted.values(model) + I(fitted.values(model)^2))

# 4. Goldfeld-Quandt test
gqtest(model, order.by = ~ x, data = adni.5)

# PCA
# Step 1: Load necessary libraries
# No extra libraries needed, but you can use ggplot2 for plotting

# Step 2: Prepare the data
# Let's assume you have a dataset 'df' with variables to run PCA on
# adni.pca <- adni.5[c('SPARE_BAG','AgeAccClockCombo','AGE','Years_to_DNAm')]
# names(adni.pca) <- c('SPARE BAG','BAGcombo','AGE','TIME')

adni.pca <- adni.5[c('SPARE_BAG','AgeAccClockCombo','AGE')]
names(adni.pca) <- c('SPARE BAG','BAGcombo','AGE')

# Step 3: Standardize the data
df_scaled <- scale(adni.pca)

# Step 4: Run PCA
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

# Step 5: Inspect results
summary(pca_result)  # To see importance of each principal component

# Eigenvectors (aka loadings)
pca_result$rotation

# Principal component scores
pca_result$x

# Step 6: Visualize the PCA
biplot(pca_result, scale = 0, xlabs=rep(".", nrow(adni.pca)))

# Step 1: Extract the proportion of variance explained by each principal component
pca_var <- pca_result$sdev^2  # Eigenvalues (variances of the principal components)
pve <- pca_var / sum(pca_var)  # Proportion of variance explained
# Optional: Add a cumulative variance line
cumulative_pve <- cumsum(pve)
lines(cumulative_pve, type = "b", col = "red", pch = 19)

# Create a dataframe with principal component index and variance explained
pca_df <- data.frame(
  PC = 1:length(pve),
  Variance_Explained = pve,
  Cumulative_Variance = cumulative_pve
)

# Plot with ggplot2
ggplot(pca_df, aes(x = PC)) +
  geom_point(aes(y = Variance_Explained), color = "blue") +
  geom_line(aes(y = Variance_Explained), color = "blue") +
  geom_point(aes(y = Cumulative_Variance), color = "red") +
  geom_line(aes(y = Cumulative_Variance), color = "red") +
  labs(title = "Scree Plot", 
       x = "Principal Component", 
       y = "Variance Explained") +
  theme_minimal()
