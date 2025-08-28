library(dplyr)

adni <- read.csv('~/Projects/hippoage/data/ADNI/20250427_adni_pheno.csv')
abc <- read.csv('~/Projects/hippoage/data/INDD/20250502_abc_pheno.csv')

boxplot(adni$M_Hippo_VOL_ASHST1,adni$M_Hippo_VOL_ASHST1.combat)

boxplot(adni$M_Hippo_VOL_ASHST1,adni$M_Hippo_VOL_ASHST1_Adj2)
boxplot(adni$M_Hippo_VOL_ASHST1_Adj1,adni$M_Hippo_VOL_ASHST1_Adj2)
boxplot(adni$M_Hippo_VOL_ASHST1_Adj1,abc$M_Hippo_VOL_ASHST1_Adj1)
t.test(adni$M_Hippo_VOL_ASHST1_Adj1,abc$M_Hippo_VOL_ASHST1_Adj1)

boxplot(adni$M_Hippo_VOL_ASHST1_Adj1,abc_w$M_Hippo_VOL_ASHST1_Adj1)
t.test(adni$M_Hippo_VOL_ASHST1_Adj1,abc_w$M_Hippo_VOL_ASHST1_Adj1)

boxplot(abc$M_Hippo_VOL_ASHST1_Adj1,abc$M_Hippo_VOL_ASHST1_Adj2) # Something is wrong with Age/ICV Adjusted

abc_w <- abc[abc$Race=="White",]

library(ggplot2)

# Assuming `adni` and `abc` are your data frames
# Combine data into a single data frame for ggplot2
df <- data.frame(
  Group = factor(c(rep("Harmonized and Adjusted ADNI", length(adni$M_Hippo_VOL_ASHST1_Adj1)),
                   rep("Adjusted ABC", length(abc$M_Hippo_VOL_ASHST1_Adj1))),
                 levels = c("Harmonized and Adjusted ADNI", "Adjusted ABC")),
  Volume = c(adni$M_Hippo_VOL_ASHST1_Adj1, abc$M_Hippo_VOL_ASHST1_Adj1)
)

# Create the ggplot boxplot
ggplot(df, aes(x = Group, y = Volume)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "Hippocampal Volume in ADNI and ABC",
    x = NULL,
    y = expression("Adjusted Hippocampal volume (mm"^3*")")
  ) +
  theme_minimal(base_size = 14)
