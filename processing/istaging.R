library(dplyr)
library(RSQLite)
library(EnvStats)
library(ggplot2)

adni <- read.csv('~/Projects/hippoage/data/ADNI/20250427_adni_pheno.csv')
abc <- read.csv('~/Projects/hippoage/data/INDD/20250530_abc_pheno.csv')

# Combine data into a single data frame for ggplot2
df <- data.frame(
  SPAREBA = c(adni$SPARE_BA, abc$SPARE_BA),
  Age = c(adni$Age_MRI, abc$Age_MRI),
  Group = factor(c(rep("ADNI", length(adni$SPARE_BA)),
                   rep("ABC", length(abc$SPARE_BA))),
                 levels = c("ADNI", "ABC"))
)

ggplot(df,aes(Age,SPAREBA)) +
  geom_point(aes(color=Group),na.rm = TRUE) + 
  labs(title="SPARE-BA vs Chronological Age", x ="Chronological Age (years)", y = "SPARE-BA (years)", fill = "Biological Age Group") +
  theme_pubclean() +
  stat_cor()

