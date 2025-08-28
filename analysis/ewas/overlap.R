library(dplyr)
library(readxl)

ewas <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/20250228_HV_EWAS.rds')
external <- read_excel('~/Downloads/41380_2019_605_MOESM4_ESM.xlsx',skip = 1,col_names = TRUE)
external_cpgs <- external$Probe_ID

# SUMMARY STATISTICS
smry <- summaryExtractTest(ewas)

# Extract probe IDs for top n CpGs
get_top_n_cgs <- function(smry, feature, top_n) {
  library(dplyr)
  pval = paste0("Pval_",feature)
  cgs <- smry %>% 
    arrange(!!sym(pval)) %>%
    slice_min(order_by = !!sym(pval), n = top_n) %>%
    select(Probe_ID)
  
  return(cgs)
}

ewas_cpgs <- get_top_n_cgs(smry,"M_Hippo_VOL_ASHST1.adjusted",60000)
overlap_cpgs <- intersect(external_cpgs,ewas_cpgs)