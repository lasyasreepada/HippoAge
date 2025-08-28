library(dplyr)

# READ DATA
ewas <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/SPAREBAG.rds')

# SUMMARY STATISTICS
smry <- summaryExtractTest(ewas)

get_fdr <- function(smry, feature) {
  library(dplyr)
  pval = paste0("Pval_",feature)
  smry$FDR <- p.adjust(smry[[pval]], method="fdr")
  return(smry)
}

smry <- get_fdr(smry, "SPAREBAG")

# Scan for top hits
print(smry %>% arrange(Pval_SPAREBAG) %>% 
  dplyr::select(Probe_ID,Est_SPAREBAG,Pval_SPAREBAG,FDR) %>% 
  filter(Pval_SPAREBAG < 0.00001),n=22)

# Extract probe IDs for top n CpGs
get_top_n_cgs <- function(smry, feature, top_n=2000) {
  library(dplyr)
  pval = paste0("Pval_",feature)
  cgs <- smry %>% 
    arrange(!!sym(pval)) %>%
    slice_min(order_by = !!sym(pval), n = top_n) %>%
    select(Probe_ID)
  
  return(cgs)
}

# Get top 2k cgs
cg20k <- get_top_n_cgs(smry,feature="M_Hippo_VOL_ASHST1_3T",top_n=20000)

# Write csv
write.csv(cg20k,'~/Projects/hippoage/data/ADNI/EWAS/20250218_HV_CGs.csv') 



  
  