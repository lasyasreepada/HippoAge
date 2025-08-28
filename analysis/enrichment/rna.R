### PWG EWAS (1A)

library(tidyverse)
library(sesame)
library(GenomicRanges)
library(SummarizedExperiment)

### Identify significant CpGs and proximal genes -- test correlation with gene expression (1B)

source("Link_Genes.R")
fdr_thresh <- .05
smry_sig <- smry %>% filter(FDR_log_pTau < fdr_thresh) %>%
  dplyr::select(Probe_ID,Est_log_pTau,FDR_log_pTau)
cg <- smry_sig$Probe_ID
bp <- 10000

cg_gene_link <- linkGenes(
  qry=cg,
  platform="EPIC",
  distance=bp
)

#read in ROSMAP methylation and gene expression data
exp_log <- readRDS("~/PART_AD/20240409_ROSMAP_GeneExpMtx_filtered_logT.rds") #expression data
se <- readRDS("~/PART_AD/ROSMAP/ROSMAP_SE.rds")
betas <- assay(se)[,colnames(exp_log)] 

cpg_gene_corr <- do.call(rbind,lapply(1:nrow(cg_gene_link),function(x) {
  cg <- cg_gene_link[x,"Probe_ID"]
  enst <- splitString(cg_gene_link[x,"ENST"],"\\.",1)
  gene <- cg_gene_link[x,"Gene"]
  if ((cg %in% rownames(betas)) & (enst %in% rownames(exp_log))) {
    res <- cor.test(betas[cg,],exp_log[enst,])
    data.frame(CpG=cg,ENST=enst,Gene=gene,Corr=res$estimate,P=res$p.value)
  } else {
    data.frame(CpG=cg,ENST=enst,Gene=gene,Corr=NA,P=NA)
  }
}))

ind <- which(!is.na(cpg_gene_corr$P)) #only the cpg-gene pairs that had non-missing probes/genes
cpg_gene_corr$FDR <- NA
cpg_gene_corr$FDR[ind] <- p.adjust(cpg_gene_corr$P[ind],method="fdr")

#Now test if these genes associate with pTau
df <- colData(se)[colnames(exp_log),]
df$log_pTau_hipp <- log(df$nft_hip + 1)

gene_nft <- do.call(rbind,lapply(1:nrow(cpg_gene_corr),function(x) {
  enst <- cpg_gene_corr[x,"ENST"]
  gene <- cpg_gene_corr[x,"Gene"]
  if(!enst %in% rownames(exp_log))  {
    return(data.frame(ENST=enst,Gene=NA,Est=NA,P=NA))
  }
  d0 <- df
  expr <- exp_log[enst,]
  d0[[enst]] <- expr
  fml <- as.formula(paste("log_pTau_hipp ~",enst))
  model <- lm(fml,data=d0)
  smry <- summary(model)$coefficients
  data.frame(
    ENST=enst,
    Gene=gene,
    Est=smry[enst,"Estimate"],
    P=smry[enst,"Pr(>|t|)"]
  )
}))
gene_nft$FDR_gene_nft <- p.adjust(gene_nft$P,method="fdr")
cpg_gene_corr <- cbind(cpg_gene_corr,gene_nft)
saveRDS(cpg_gene_corr,file="~/PART_AD/PWG/PWG_PART_pTau_EWAS_Hits_CpG_Gene_Corr.rds")
```

### Test correlation of EWAS hits with all genes (1C,S2B)
betas_sig <- betas[intersect(cg,rownames(betas)),] #betas subsetted for the EWAS hits

idx = seq_len(nrow(betas_sig))
idx_exp <- seq_len(nrow(exp_log))
cpg_gene_corr_all <- do.call(rbind,lapply(idx,function(x) tibble(i=rep(x,nrow(exp_log)),j=idx_exp))) #indices for each cg in betas mtx, with each gene in expression matrix

tmp <- do.call(rbind,mclapply(seq_len(nrow(cpg_gene_corr_all)), function(m) {
  beta_vals <- betas_sig[cpg_gene_corr_all[m,]$i,] #get betas for the cg at row m
  expr_vals <- exp_log[cpg_gene_corr_all[m,]$j,]
  res <- cor.test(beta_vals,expr_vals, use="na.or.complete")
  data.frame(Cor=res$estimate,P=res$p.value)
},mc.cores=40)) #cpg - gene correlations for all
tmp$FDR <- p.adjust(tmp$P,method="fdr")

cpg_gene_corr_all <- cbind(cpg_gene_corr_all,tmp)
cpg_gene_corr_all$CpG <- rownames(betas)[cpg_gene_corr_all$i]
cpg_gene_corr_all$ENST <- rownames(exp_log)[cpg_gene_corr_all$j]

#match gene names to transcript ids
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

gene_names <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                    filters = "ensembl_gene_id",
                    values = rownames(exp_log),
                    mart = ensembl.con)

cpg_gene_corr_all$Gene <- gene_names$external_gene_name[match(cpg_gene_corr_all$ENS,gene_names$ensembl_gene_id)] 
cpg_gene_corr_all <- cpg_gene_corr_all[,c("CpG","ENST","Gene","Cor","P","FDR")]
saveRDS(cpg_gene_corr_all,file="~/PART_AD/PWG/PWG_PART_pTau_EWAS_Hits_CpG_Gene_Corr_allGenes.rds")
