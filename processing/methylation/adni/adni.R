library(dplyr)
library(sesame)
library(snow)
library(stats)
library(RSQLite)
library(BiocParallel)
library(SummarizedExperiment)
library(lubridate)

# DATA IO
idat_dir_adni <- "/Users/sreepada/Library/CloudStorage/Box-Box/DataRepository/ADNI/Methylation/ADNI_iDAT_files/"

# ADNI
adni_betas <- openSesame(idat_dir_adni, func = getBetas)
adni_qcs <- openSesame(idat_dir_adni, prep="", func=sesameQC_calcStats, funs="detection")

# Get names
adni_samples <- colnames(adni_betas)

# ADNI
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in subject information (coldata)
adni_m_annot_file <- "data/ADNI/methylation/ADNI_DNA_Methylation_SampleAnnotation_20170530.xlsx"
adni_m_annot <- readxl::read_xlsx(adni_m_annot_file) # ADNI

# Put DNAm assay date in datetime format
adni_m_annot$Edate <- as.Date(adni_m_annot$Edate, format = "%m/%d/%Y")

# Read in Demographics Data
adniDems <- read.csv("~/Projects/havana/data/ADNI/Subjects/ADNI_Demographics.csv")
adniDems$PTDOBDD <- 1
adniDems$PTDOB <- as.Date(with(adniDems, paste(PTDOBMM, PTDOBDD, PTDOBYY,sep="/")), format="%m/%d/%Y")
adniDems <- subset(adniDems, select = c(RID,PTDOB,PTGENDER))

adniDems <- adniDems %>% 
  group_by(RID) %>%
  arrange(PTDOB) %>%
  filter(row_number()==1)

# Join demographics with annotation table
adni_m_annot <- left_join(adni_m_annot,adniDems,by="RID")

# Compute chronological Age of DNAm assay
adni_m_annot$Age <- time_length(adni_m_annot$Edate - adni_m_annot$PTDOB,"years")

# USE ONLY TO UPDATE COLDATA IN AN EXISTING RDS/SE
# adniDNAm <- readRDS("~/Projects/havana/data/ADNI/Methylation/methylation.rds")
# adniDNAmInfo <- as.data.frame(colData(adniDNAm))
# adni_betas <- assays(adniDNAm)[[1]]
# adni_qcs <- metadata(adniDNAm)
# adni_samples <- as.data.frame(colnames(adni_betas))
colnames(adni_samples) <- "barcodes"

# Remove participants missing age
missingAge <- adni_samples[!adni_samples$barcodes %in% adni_m_annot$barcodes,]
adni_samples <- adni_samples[!adni_samples$barcodes %in% missingAge,]
adni_samples <- as.data.frame(adni_samples)
colnames(adni_samples) <- "barcodes"
betas <- adni_betas[ , !colnames(adni_betas) %in% missingAge]

# Prep DNAm samples to be merged with annotation
names(adni_m_annot)[names(adni_m_annot) == "PTGENDER"] <- "Sex"
adni_m_annot$Sex <- ifelse(adni_m_annot$Sex == 1,"Male","Female")

# Match methylation samples with annotation based on 'barcodes'
adni_colData <- as.data.frame(full_join(adni_samples, adni_m_annot, by="barcodes"))
names(adni_colData)[names(adni_colData) == "barcodes"] <- "Sample"

# Save as SummarizedExperiment
se <- SummarizedExperiment(assays = betas, colData=adni_colData, metadata = adni_qcs)

# Checkpoint
identical(colnames(assay(se)), rownames(colData(se)))

# Save to RDS file
saveRDS(se, "data/ADNI/methylation/methylation.rds")

# Useful vignettes and checkpoints
# https://www.bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/inferences.html


