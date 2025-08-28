library(dplyr)
library(sesame)
library(knowYourCG)
library(SummarizedExperiment)
library(GenomicRanges)
library(devtools)

setwd('/Users/lasyasreepada/projects/hippoage')
source("analysis/enrichment/linkgenes.R")

# Read EN feature importances
spareba <- read.csv('~/Projects/hippoage/data/model_coefs_SPARE_BA.csv')
spareba <- spareba[(!spareba$Feature.Names == 'Intercept' & !spareba$Coefficients==0),]

hyper <- spareba[spareba$Coefficients > 0,]
hypo <- spareba[spareba$Coefficients < 0,]

# Get overlap
cpgs <- spareba$Feature.Names

# KYCG ENRICHMENT

# Link genes
bp <- 10000
cg_gene_link_sba <- linkGenes(
  qry=cpgs,
  platform="EPIC",
  distance=bp
)

devtools::install_github("zhou-lab/knowYourCG")
kycg <- testEnrichment(cpgs, platform = 'EPIC')

# GENE ONTOLOGY
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
write.csv(enrich_df,'~/Projects/hippoage/data/ADNI/SPAREBA_GO.csv')

# To Excel
printEnrich(enriched, outFile = "excel")

