library(dplyr)
library(sesame)
library(knowYourCG)
library(SummarizedExperiment)
library(GenomicRanges)

setwd('/Users/lasyasreepada/projects/hippoage')
source("analysis/enrichment/linkgenes.R")

# Read EWAS
ewas <- readRDS('data/ADNI/EWAS/20250308_HV_ADNI_FINAL.rds') 
smry <- summaryExtractTest(ewas)
smry$FDR_M_Hippo_VOL_ASHST1.combat <- p.adjust(smry$Pval_M_Hippo_VOL_ASHST1.combat,method="fdr")

threshold <- 5e-4

smry_sig <- smry %>% filter(Pval_M_Hippo_VOL_ASHST1.combat < threshold) %>%
  dplyr::select(Probe_ID,Est_M_Hippo_VOL_ASHST1.combat,Pval_M_Hippo_VOL_ASHST1.combat,FDR_M_Hippo_VOL_ASHST1.combat)

smry %>% arrange(Pval_M_Hippo_VOL_ASHST1.combat) %>% filter(Pval_M_Hippo_VOL_ASHST1.combat < threshold) %>%
  dplyr::select(Probe_ID,Est_M_Hippo_VOL_ASHST1.combat,Pval_M_Hippo_VOL_ASHST1.combat,FDR_M_Hippo_VOL_ASHST1.combat)

# Read EN feature importances
hv <- read.csv('~/Projects/hippoage/data/ADNI/20250305_HV_EN_Features.csv')
kycg <- annoProbes(hv$name,platform = 'EPIC')

# Get overlap
overlap <- intersect(smry_sig$Probe_ID,hv$name)
overlap <- hv$name
  
# Link genes
bp <- 10000
cg_gene_link <- linkGenes(
  qry=overlap,
  platform="EPIC",
  distance=bp
)

library(enrichR)

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023")

input <- cg_gene_link$Gene

if (websiteLive) {
  enriched <- enrichr(input, dbs)
}

biol <- enriched[["GO_Biological_Process_2023"]]
mol <- enriched[["GO_Molecular_Function_2023"]]
cell <- enriched[["GO_Cellular_Component_2023"]]

if (websiteLive) {
  plotEnrich(biol, showTerms = 40, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(mol, showTerms = 8, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(cell, showTerms = 8, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

# CSV
enriched_df <- rbind(biol,mol, cell)
write.csv(enriched_df,'~/Projects/hippoage/data/ADNI/20250305_HV_Enriched.csv',row.names = FALSE)
