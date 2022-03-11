#Jenny Smith
#Feb 22, 2022

library(dplyr)
library(readr)
library(stringr)
library(magrittr)


#### Request

## Read in the counts data 
genome <- "GRCh38"
## GRCh38 
current_files <- dir(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/"))
# current_files


if(genome=="GRCh38"){
  grch38_cts_file <- grep("_RBD_.+scaledTPM_counts.RDS", current_files, value=TRUE)
  cts_grch38 <- readRDS(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/",grch38_cts_file))
  
  ### TPM
  grch38_TPM_file <- grep("_RBD_.+Abundance_TPM", current_files, value=TRUE)
  TPM_grch38 <- readRDS(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/",grch38_TPM_file))
}

# head(cts_grch38[,1:5])
# head(TPM_grch38[,1:5])


## Sample manifest and clinical data
sample_info <- read.csv(file.path(TARGET, "SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_08.12.21.csv")) 

dim(sample_info)


# merged <- read.csv(file.path(CDE, "Merged/TARGET_AML_0531_1031_merged_CDEs_05.21.21.csv"), 
#                    na.strings = c("N/A","#N/A","NA","^$", "^\\.$")) 
# 
# inelig <- merged %>% 
#   filter(Eligibility_Comments == "remove") %>% 
#   pull(USI)


## Subset the Sample Manifest
ds_normal_samples <- sample_info %>% 
  filter(grepl("DS|TMD|^NBM|CD34_PB", Group)) %>% 
  filter(Sample %in% colnames(cts_grch38))


dim(ds_normal_samples)
table(ds_normal_samples$Group)
table(ds_normal_samples$Time_point)


#Save the Data 
columns <- c("gene_id","gene_name",ds_normal_samples$Sample)
# write.csv(ds_normal_samples,
#           file.path(SCRATCH,"jlsmith3/TARGET_DS-AML_TMD_HealthyControls_GRCh38_GencodeV29_RNA-seq_Manifest.csv"), row.names = FALSE)
# write.csv(cts_grch38[,columns],
#           file.path(SCRATCH,"jlsmith3/TARGET_DS-AML_TMD_HealthyControls_GRCh38_GencodeV29_dupGenesRemoved_FractionalCounts.csv"), row.names = FALSE)
# write.csv(TPM_grch38[,columns],
#           file.path(SCRATCH,"jlsmith3/TARGET_DS-AML_TMD_HealthyControls_GRCh38_GencodeV29_dupGenesRemoved_TPM.csv"), row.names = FALSE)
# 

