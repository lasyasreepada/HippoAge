# Load libraries
library(hgu219.db)
library(biomaRt)
library(dplyr)
library(AnnotationDbi) 
library(SummarizedExperiment)
library(sesame)

sesameDataCache()

source("~/Projects/hippoage/analysis/enrichment/linkgenes.R")

# Get all probe IDs
all_probes <- keys(hgu219.db, keytype = "PROBEID")

# Map to gene symbols and Entrez IDs
probe2gene <- AnnotationDbi::select(hgu219.db,
                     keys = all_probes,
                     columns = c("SYMBOL", "ENTREZID"),
                     keytype = "PROBEID")

# Step 2: Convert ENTREZ to character if needed
probe2gene$ENTREZID <- as.character(probe2gene$ENTREZID)

# Step 3: Query biomaRt to map ENTREZ to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_mapping <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "entrezgene_id",
  values = probe2gene$ENTREZID,
  mart = ensembl
)

# Step 4: Make ENTREZID a character for joining
ensembl_mapping$entrezgene_id <- as.character(ensembl_mapping$entrezgene_id)

# Step 5: Join your data to Ensembl mapping
final <- left_join(probe2gene, ensembl_mapping,
                   by = c("ENTREZID" = "entrezgene_id"), 
                   relationship = 'many-to-many') 
final <- final %>% distinct()

# Read ADNI Gene Expression Data
genex <- read.csv('~/Projects/hippoage/data/ADNI/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv',header = FALSE)

# Phenotypes
pheno <- genex[1:8,]
pheno <- dplyr::select(pheno, -c('V2', 'V3','V748'))
pheno <- setNames(data.frame(t(pheno[,-1])), pheno[,1])
rownames(pheno) <- pheno$SubjectID

# Gene expression values
exp <- genex[10:nrow(genex),-ncol(genex)]
rownames(exp) <- exp[[1]]
exp <- exp[, -1]
exp <- exp[,3:ncol(exp)]
colnames(exp) <- pheno$SubjectID

# Create Summarized Experiment
# se <- SummarizedExperiment(assays = exp, colData=pheno)

# link keys to ensembl IDs
map <- genex[9:nrow(genex),c(1:3,ncol(genex))]
colnames(map) <- c(as.character(map[1, -ncol(map)]),'Description')
map <- map[-1, ]
rownames(map) <- NULL
map$LocusLink <- sub("^LOC", "", map$LocusLink)
map$LocusLink[map$LocusLink == ""] <- NA
map <- map %>%
  dplyr::rename(PROBEID = ProbeSet,
                ENTREZID = LocusLink,
                SYMBOL = Symbol)

map <- left_join(map, final[,c('PROBEID','ensembl_gene_id')], by='PROBEID')

# Read EN feature importances
adni <- read.csv('~/Projects/hippoage/data/model_coefs_M_Hippo_VOL_ASHST1_Adj1.csv')
abc <- read.csv('~/Projects/hippoage/data/model_coefs_abc_M_Hippo_VOL_ASHST1_Adj1.csv')

adni <- adni[(!adni$Feature.Names == 'Intercept' & !adni$Coefficients==0),]
abc <- abc[(!abc$Feature.Names == 'Intercept' & !abc$Coefficients==0),]

# Get overlap
cpgs <- intersect(adni$Feature.Names,abc$Feature.Names)

# KYCG ENRICHMENT

# Link genes
bp <- 10000
cg_gene_link <- linkGenes(
  qry=cpgs,
  platform="EPIC",
  distance=bp
)

cg_gene_link_adni <- linkGenes(
  qry=adni$Feature.Names,
  platform="EPIC",
  distance=bp
)

cg_gene_link_abc <- linkGenes(
  qry=abc$Feature.Names,
  platform="EPIC",
  distance=bp
)

overlap_genes <- intersect(cg_gene_link_adni$Gene,cg_gene_link_abc$Gene)
overlap_df_adni <- cg_gene_link_adni[cg_gene_link_adni$Gene %in% overlap_genes,]
overlap_df_abc <- cg_gene_link_abc[cg_gene_link_abc$Gene %in% overlap_genes,c('Probe_ID','Gene')]
overlap_df <- merge(overlap_df_adni,overlap_df_abc,by=c('Gene'))

# GENE ONTOLOGY
library(enrichR)
websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023")

input <- 'TREML2'

if (websiteLive) {
  enriched <- enrichr(input, dbs)
}

biol <- enriched[["GO_Biological_Process_2023"]]
mol <- enriched[["GO_Molecular_Function_2023"]]
cell <- enriched[["GO_Cellular_Component_2023"]]

if (websiteLive) {
  plotEnrich(biol, showTerms = 25, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(mol, showTerms = 5, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(cell, showTerms = 3, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

# Write CSV
enrich_df <- rbind(biol,mol, cell)
write.csv(enrich_df,'~/Projects/hippoage/data/HV_GO.csv')

# Remove versioning from Ensembl IDs
splitString <- function (x, s, i) {
  l <- strsplit(x, s)
  vapply(l, function(x) x[[i]], character(1))
}

# CORRELATION WITH BLOOD GENE EXPRESSION
cg_gene_link$ENST <- sub("\\..*", "", cg_gene_link$ENST)

exp$PROBEID <- rownames(exp)
exp_filtered <- left_join(map,exp,by='PROBEID')
exp_filtered <- exp_filtered %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

rownames(exp_filtered) <- exp_filtered$ensembl_gene_id

# Read SE
se <- readRDS('~/Projects/hippoage/data/ADNI/20250612_ADNI_SE.rds')
se <- se[,se$PTID %in% colnames(exp_filtered)]
betas <- assay(se)
colnames(betas) <- se$PTID
pheno <- as.data.frame(colData(se))

# Set matrices
meth_matrix <- betas
expr_matrix <- exp_filtered[,colnames(meth_matrix)]
expr_matrix_numeric <- as.data.frame(sapply(expr_matrix, as.numeric))
expr_matrix_numeric <- log(expr_matrix_numeric)
rownames(expr_matrix_numeric) <- rownames(expr_matrix)
colnames(expr_matrix_numeric) <- colnames(expr_matrix)

# meth_matrix_numeric <- na.omit(meth_matrix)
# expr_matrix_numeric <- na.omit(expr_matrix_numeric)


# Manual Correlation

# BCL11B
cg1a <- 'cg23479730'
cg1b <- 'cg05100282'
enst1a <- cg_gene_link[cg_gene_link$Probe_ID==cg1a,'ENST']
enst1b <- enst1a

# TREML2
cg2 <- 'cg03526776'
enst2 <- cg_gene_link[cg_gene_link$Probe_ID==cg2,'ENST']

# Corr
cor.test(as.numeric(meth_matrix[cg1a,]),as.numeric(expr_matrix_numeric[enst1a,]))
cor.test(as.numeric(meth_matrix[cg1b,]),as.numeric(expr_matrix_numeric[enst1b,]))
cor.test(as.numeric(meth_matrix[cg2,]),as.numeric(expr_matrix_numeric[enst2,]))

# FUNCTION TO COMPUTE GENE OVERLAP LIKELIHOOD
# This function calculates the probability of observing the observed overlap
# or more extreme overlap between two gene sets by chance

compute_gene_overlap_likelihood <- function(set1_genes, set2_genes, background_genes = NULL) {
  # If background is not provided, use union of both sets plus some common background
  if (is.null(background_genes)) {
    background_genes <- unique(c(set1_genes, set2_genes))
  }
  
  # Remove any NA or empty values
  set1_genes <- set1_genes[!is.na(set1_genes) & set1_genes != ""]
  set2_genes <- set2_genes[!is.na(set2_genes) & set2_genes != ""]
  background_genes <- background_genes[!is.na(background_genes) & background_genes != ""]
  
  # Calculate overlap
  overlap_genes <- intersect(set1_genes, set2_genes)
  n_overlap <- length(overlap_genes)
  
  # Set sizes
  n_set1 <- length(set1_genes)
  n_set2 <- length(set2_genes)
  n_background <- length(background_genes)
  
  # Perform hypergeometric test
  # phyper(q, m, n, k, lower.tail = FALSE) where:
  # q = number of white balls drawn (overlap - 1, for P(X >= overlap))
  # m = number of white balls in urn (set1 size)
  # n = number of black balls in urn (background - set1 size)
  # k = number of balls drawn (set2 size)
  
  p_value <- phyper(n_overlap - 1, n_set1, n_background - n_set1, n_set2, lower.tail = FALSE)
  
  # Calculate enrichment ratio (observed/expected)
  expected_overlap <- (n_set1 * n_set2) / n_background
  enrichment_ratio <- n_overlap / expected_overlap
  
  # Create result object
  result <- list(
    set1_size = n_set1,
    set2_size = n_set2,
    background_size = n_background,
    overlap_size = n_overlap,
    overlap_genes = overlap_genes,
    expected_overlap = expected_overlap,
    enrichment_ratio = enrichment_ratio,
    p_value = p_value,
    significant = p_value < 0.05
  )
  
  return(result)
}

# Example usage:
result <- compute_gene_overlap_likelihood(
  set1_genes = unique(cg_gene_link_adni$Gene),
  set2_genes = unique(cg_gene_link_abc$Gene),
  background_genes = union(cg_gene_link_adni$Gene, cg_gene_link_abc$Gene)
)

print(paste("Overlap size:", result$overlap_size))
print(paste("Expected overlap:", round(result$expected_overlap, 2)))
print(paste("Enrichment ratio:", round(result$enrichment_ratio, 2)))
print(paste("P-value:", format(result$p_value, scientific = TRUE)))
print(paste("Significant:", result$significant))
