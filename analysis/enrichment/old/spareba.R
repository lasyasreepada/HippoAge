library(dplyr)
library(sesame)
library(knowYourCG)
library(SummarizedExperiment)
library(GenomicRanges)

setwd('/Users/lasyasreepada/projects/hippoage')
source("analysis/enrichment/linkgenes.R")

# Read EWAS
smry <- readRDS('data/ADNI/EWAS/20250305_SPAREBA_EWAS2.rds') 
smry$FDR_SPAREBA <- p.adjust(smry$Pval_SPAREBA,method="fdr")

threshold <- 5e-4

smry_sig <- smry %>% filter(Pval_SPAREBA < threshold) %>%
  dplyr::select(Probe_ID,Est_SPAREBA,Pval_SPAREBA,FDR_SPAREBA)

smry %>% arrange(Pval_SPAREBA) %>% filter(Pval_SPAREBA < threshold) %>%
  dplyr::select(Probe_ID,Est_SPAREBA,Pval_SPAREBA,FDR_SPAREBA)

# Read EN feature importances
spareba <- read.csv('~/Projects/hippoage/data/ADNI/20250305_SPAREBA_EN_Features.csv')

# Get overlap
overlap <- union(smry_sig$Probe_ID,spareba$name)
overlap <- spareba$name

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
  plotEnrich(biol, showTerms = 20, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(mol, showTerms = 8, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(cell, showTerms = 5, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

# Write CSV
enrich_df <- rbind(biol,mol, cell)
write.csv(enrich_df,'~/Projects/hippoage/data/ADNI/20250305_SPAREBA_Enriched.csv')

# To Excel
printEnrich(enriched, outFile = "excel")

