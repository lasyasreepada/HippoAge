# CHARACTERISTICS
myVars <- c('Age','Sex','apoe','Demographic.Years.of.Education','Neuropsych.CDR.SOB')
catVars <- c('Sex','apoe','Demographic.Years.of.Education')

# FACTOR
aibl$Diagnosis <- factor(aibl$Diagnosis, levels = c('HC','MCI','AD')) 
aibl$Demographic.Years.of.Education <- factor(aibl$Demographic.Years.of.Education, levels = c('0-6','7-8','9-12','13-15','15+')) 

# TABLE 1
tab1_dx <- CreateTableOne(vars = myVars, strata = 'Diagnosis', data = aibl, factorVars = catVars)

# GET UNIQUE IDS
ids_DNA <- as.integer(unique(pheno$AIBL_ID))
ids_mri <- as.integer(unique(aibl_mri$Subject))

# INTERSECT
length(intersect(ids_DNA,ids_mri))

# VENN
venn_data <- list(MRI = ids_mri, DNAM = ids_DNA)
list <- setdiff(ids_DNA,intersect(ids_mri,ids_DNA))

# Plot the Venn diagram
ggVennDiagram(venn_data) +
  scale_fill_gradient(
    low = "lightblue",  # Color for low values
    high = "darkblue",  # Color for high values
  ) +
  ggtitle('Venn Diagram of Imaging and Methylation overlap')
