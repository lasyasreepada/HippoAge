library(dplyr)
library(ggplot2)
library(forestplot)
library(boot)

spareba <- read.csv('~/Projects/hippoage/data/ADNI/20250305_SPAREBA_EN_Features.csv')
hv <- read.csv('~/Projects/hippoage/data/ADNI/20250305_HV_EN_Features.csv')

spareba <- as.data.frame(spareba)
hv <- as.data.frame(hv)

names(spareba) <- c('feature','importance')
names(hv) <- c('feature','importance')

# Prepare the data for plotting
spareba_plot <- spareba %>%
  mutate(feature = factor(feature, levels = feature[order(importance)]))

# Prepare the data for plotting
hv_plot <- hv %>%
  mutate(feature = factor(feature, levels = feature[order(importance)]))


# Create the plot
ggplot(spareba_plot, aes(x = importance, y = feature)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.3f", importance)), hjust = -0.2, size = 3.5) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Feature Importances for SPARE Brain Age Clock",
    x = "Importance"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)))

# Create the plot
ggplot(hv_plot, aes(x = importance, y = feature)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.3f", importance)), hjust = -0.2, size = 3.5) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = "Feature Importances for Hippocampal Volume Clock",
    x = "Importance"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)))




