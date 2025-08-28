library(dplyr)
library(sesame)
library(knowYourCG)
library(SummarizedExperiment)
library(GenomicRanges)
library(devtools)

setwd('/Users/lasyasreepada/projects/hippoage')
source("analysis/enrichment/linkgenes.R")

# Read EN feature importances
adni <- read.csv('~/Projects/hippoage/data/model_coefs_M_Hippo_VOL_ASHST1_Adj1.csv')
abc <- read.csv('~/Projects/hippoage/data/model_coefs_abc_M_Hippo_VOL_ASHST1_Adj1.csv')

adni <- adni[(!adni$Feature.Names == 'Intercept' & !adni$Coefficients==0),]
abc <- abc[(!abc$Feature.Names == 'Intercept' & !abc$Coefficients==0),]

adni_hyper <- adni[adni$Coefficients > 0,'Feature.Names']
adni_hypo <- adni[adni$Coefficients < 0,'Feature.Names']

abc_hyper <- abc[abc$Coefficients > 0,'Feature.Names']
abc_hypo <- abc[abc$Coefficients < 0,'Feature.Names']

adni$Coefficients <- abs(adni$Coefficients)
abc$Coefficients <- abs(abc$Coefficients)

# top 25%
top_n_adni <- adni %>%
  slice_max(order_by = Coefficients, prop = 0.27)

top_n_abc <- abc %>%
  slice_max(order_by = Coefficients, prop = 0.27)

# Get overlap
cpg_overlap <- intersect(adni$Feature.Names, abc$Feature.Names)

# KYCG ENRICHMENT

# Link genes
bp <- 10000
# top_n_adni_link <- linkGenes(
#   qry=top_n_adni$Feature.Names,
#   platform="EPIC",
#   distance=bp
# )
# 
# top_n_abc_link <- linkGenes(
#   qry=top_n_abc$Feature.Names,
#   platform="EPIC",
#   distance=bp
# )

adni_link <- linkGenes(
  qry=adni$Feature.Names,
  platform="EPIC",
  distance=bp
)

abc_link <- linkGenes(
  qry=abc$Feature.Names,
  platform="EPIC",
  distance=bp
)

overlap <- intersect(adni_link$Gene,abc_link$Gene)

devtools::install_github("zhou-lab/knowYourCG")
kycg <- testEnrichment(adni_link$Feature.Names, platform = 'EPIC')

# GENE ONTOLOGY
library(enrichR)
websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023")

input <- overlap

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
enrich_overlap <- rbind(biol,mol, cell)
write.csv(enrich_df,'~/Projects/hippoage/data/INDD/HV1_GO.csv')

# To Excel
# printEnrich(enriched, outFile = "excel")

