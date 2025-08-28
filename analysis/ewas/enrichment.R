library(dplyr)
library(sesame)

# READ DATA
ewas <- readRDS('~/Projects/hippoage/data/ADNI/EWAS/SPAREBA.rds')

# Extract summary
smry <- summaryExtractTest(ewas)
smry$FDR <- p.adjust(smry$Pval_SPAREBA, method="fdr")
smry %>% arrange(Pval_SPAREBA) %>% dplyr::select(Probe_ID,Est_SPAREBA,Pval_SPAREBA,FDR) 

# PLOT RESULTS

# Volcano
smry$estimate <- smry$Est_SPAREBA
p1 <- KYCG_plotVolcano(smry,label_by = "Probe_ID",alpha = 0.2)

pdf("~/Projects/hippoage/illustrations/20250228_SPAREBAG_volcano.pdf")
print(p1)
dev.off()

# Manhattan
age_p <- setNames(-log10(smry$Pval_SPAREBA),smry$Probe_ID)
p2 <- KYCG_plotManhattan(age_p,platform = "EPIC",title = 'Manhattan Plot of EWAS on SPAREBA',label_min=5)
pdf("~/Projects/hippoage/illustrations/20250301_SPAREBA_Manhattan.pdf",width = 11, height=8.5)
print(p2)
dev.off()

# TEST ENRICHMENT
qry_hyper <- smry %>% 
  filter(Pval_SPAREBA < .00001, Est_SPAREBA > 0) %>% 
  pull(Probe_ID)

res_hyper <- testEnrichment(qry_hyper,platform = "EPIC")
res_hyper$log10.p.value.neg <- -(res_hyper$log10.p.value)

p3 <- KYCG_plotEnrichAll(res_hyper)
pdf("~/Projects/hippoage/illustrations/20250228_SPAREBAG_enrichment.pdf")
print(p3)
dev.off()

p4 <- KYCG_plotDot(res_hyper,y="log10.p.value.neg",order_by = 'p.value')
pdf("~/Projects/hippoage/illustrations/20250228_SPAREBAG_dot.pdf")
print(p4)
dev.off()

