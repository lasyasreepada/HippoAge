library(dplyr)

rosmap <- read.csv('~/Projects/havana/data/ROSMAP/Methylation/DNAmClocks.csv')
ids <- unique(rosmap$subject_id)

write.csv(ids,'~/Projects/havana/data/ROSMAP/Methylation/ROSMAP_DNA_IDs.csv')