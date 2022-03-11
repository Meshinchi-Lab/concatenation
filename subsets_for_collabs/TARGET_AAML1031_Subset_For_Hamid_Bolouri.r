#Jenny Smith
#1/18/20
#Purpose: Provide counts data to Hamid Bolouri

setwd(file.path(SCRATCH,"jlsmith3"))

library(dplyr)
library(tibble)

merged <- read.csv(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_12.09.20.csv"))
  
merged <- merged %>% 
  dplyr::filter(!is.na(Reg.),
                !is.na(USI),
                USI != "Unknown") %>% 
  filter(Eligibility_Comments!="remove") %>%
  dplyr::select(-Group) %>% 
  mutate(Reg.=as.character(Reg.))

dim(merged) #2217  150

manifest <- read.csv(file.path(TARGET, "SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_10.08.20.csv"))

# head(manifest)
dim(manifest)
# intersect(colnames(manifest), colnames(merged))



######### 

manifest.forHamid <- manifest %>% 
  filter(grepl("NBM|CD34_PB|diagnostic|relapse", Time_point), 
         grepl("AML|NBM|CD34_PB", Group)) %>% 
  dplyr::select(-one_of(c("Primary.Fusion", "Primary.CNV", "Additional.Fusions.CNV",
                   "PATIENT_ID_Original", "Final_Patient_ID"))) 

dim(manifest.forHamid) #2018   10
table(manifest.forHamid$Group)
table(manifest.forHamid$Time_point)
# write.csv(manifest.forHamid, "TARGET_AML_RBD_Diagnostic_Relapse_RNAseq_Samples.csv", row.names = FALSE)

temp <- manifest.forHamid %>% 
  dplyr::filter(AML_Subtype!="NBM", AML_Subtype != "CD34_PB") %>% 
  dplyr::select(USI,Sample,AML_Subtype) 

# dplyr::filter(temp, !USI %in% merged$USI)
  # left_join(., dplyr::select(merged, USI, Eligibility_Comments, Eligibility), 
  #           by="USI") %>% 
  # filter(Eligibility_Comments=="remove")

# write.csv(temp,"TARGET_AML_Ineligables_01.22.20.csv", row.names = F)
# View(temp)


CDE.forHamid <- merged %>% 
  filter(USI %in% manifest.forHamid$USI) %>% 
  select(-Reg.)

dim(CDE.forHamid) #1517  149
table(CDE.forHamid$Protocol)
# write.csv(CDE.forHamid, "TARGET_AML_RBD_Diagnostic_Relapse_RNAseq_Samples_CDE.csv", row.names = FALSE)


######## Counts Data

#BCCA
rbd_counts <- readRDS(file.path(PROJHOME,"0000.00.03_ExpressionMatrices/TARGET_AML_MPN_DS_NBM_2646Samples_Ribodepleted_RNAseq_geneLevel_dupGenesRemoved_FractionalCounts.RDS")) %>% 
  column_to_rownames("Gene")

dim(rbd_counts) #51573  2418
head(rbd_counts[,1:5])

in_samps <- intersect(colnames(rbd_counts), manifest.forHamid$Sample)

write.csv(rbd_counts[,in_samps], "TARGET_AML_NBM_Ribodepleted_RNAseq_geneLevel_dupGenesRemoved_BCCA_FractionalCounts.csv")


rbd_TPM <- readRDS(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/TARGET_AML_MPN_DS_NBM_2646Samples_Ribodepleted_RNAseq_geneLevel_dupGenesRemoved_TPM.RDS")) %>% 
  column_to_rownames("Gene")

dim(rbd_TPM)
head(rbd_TPM[,1:5])

# write.csv(rbd_TPM[,in_samps], "TARGET_AML_NBM_Ribodepleted_RNAseq_geneLevel_dupGenesRemoved_BCCA_TPM.csv")



#Kallisto 


kallisto_gene.cts <- readRDS(file.path(PROJHOME,"0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/TARGET_AML_RBD_Dx_Rlps_NBM_MPN_Kallisto_Quant_GeneLevel_scaledTPM_counts.RDS"))

dim(kallisto_gene.cts)
head(kallisto_gene.cts[,1:5])

kall_in_samps <- intersect(colnames(kallisto_gene.cts), manifest.forHamid$Sample)

# write.csv(kallisto_gene.cts[,kall_in_samps], "TARGET_AML_NBM_Kallisto_Quant_GeneLevel_scaledTPM_counts.csv")


kallisto_gene.TPM <- readRDS(file.path(PROJHOME,"0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/TARGET_AML_RBD_Dx_Rlps_NBM_MPN_Kallisto_Quant_GeneLevel_Abundance_TPM.RDS"))

dim(kallisto_gene.TPM)
head(kallisto_gene.TPM[,1:5])

# write.csv(kallisto_gene.TPM[,kall_in_samps], "TARGET_AML_NBM_Kallisto_Quant_GeneLevel_Abundance_TPM.csv")


kallisto_tx.cts <- readRDS(file.path(PROJHOME,"0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/TARGET_AML_RBD_Dx_Rlps_NBM_MPN_Kallisto_Quant_TranscriptLevel_scaledTPM_counts.RDS"))

dim(kallisto_tx.cts)
head(kallisto_tx.cts[,1:5])

# write.csv(kallisto_tx.cts[,kall_in_samps], "TARGET_AML_NBM_MPN_Kallisto_Quant_TranscriptLevel_scaledTPM_counts.csv")



kallisto_tx.TPM <- readRDS(file.path(PROJHOME,"0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/TARGET_AML_RBD_Dx_Rlps_NBM_MPN_Kallisto_Quant_TranscriptLevel_Abundance_TPM.RDS"))

dim(kallisto_tx.TPM)
head(kallisto_tx.TPM[,1:5])

# write.csv(kallisto_tx.TPM[,kall_in_samps], "TARGET_AML_NBM_Kallisto_Quant_TranscriptLevel_Abundance_TPM.csv")









