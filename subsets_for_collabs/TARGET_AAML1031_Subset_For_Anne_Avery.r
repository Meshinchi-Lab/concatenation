#Jenny Smith
#Dec 4, 2019

library(dplyr)
library(readr)
library(stringr)
library(magrittr)


#### Request

##  Dec 3, 2021
# Hi Jenny,
# I have a quick question and Tim said you'd be the one to ask.  
# We're working on our FUS-ERG paper at the moment and we've found the kids present with an immune evasion phenotype 
# that mirrors that of adult AML post HSCT relapse.  Is there transcript data on any TARGET samples that have relapsed post-transplant? 
# Additionally, are there any other timepoints for the FUS-ERG patients besides diagnosis where RNAseq was run?  Thanks for any help you can provide!


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
diagnostic_normal_samples <- sample_info %>% 
  filter(grepl("^AML$|^Flow|^NBM", Group), 
         grepl("diagnostic|^NBM", Time_point)) %>% 
  mutate_at(vars(Group), ~case_when(
    .=="FlowSorted" ~ "AML", 
    TRUE ~ .)) %>% 
  filter(!grepl("replicate", Sample)) %>% 
  filter(Sample %in% colnames(cts_grch38))




dim(diagnostic_normal_samples)
table(diagnostic_normal_samples$Group)
table(diagnostic_normal_samples$Time_point)


#Save the Data 
columns <- c("gene_id","gene_name",diagnostic_normal_samples$Sample)
# write.csv(diagnostic_normal_samples,
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_GRCh38_GencodeV29_RNA-seq_Manifest.csv"), row.names = FALSE)
# write.csv(cts_grch38[,columns],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_GRCh38_GencodeV29_dupGenesRemoved_FractionalCounts.csv"), row.names = FALSE)
# write.csv(TPM_grch38[,columns],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_GRCh38_GencodeV29_dupGenesRemoved_TPM.csv"), row.names = FALSE)

cts_grch38[,columns[1:10]] %>%  head()
cts_grch38[,columns] %>%  dim()

#https://stackoverflow.com/questions/55787464/im-having-trouble-loading-my-csv-from-dropbox-to-rstudio-cloud-can-you-please
# url <- "https://www.dropbox.com/s/0ionb8ycgbxqed9/CSU_Canine_AML_RNAseq_STAR_ReadsPerGene_counts.csv?dl=1"
# canine_cts <- data.table::fread(url, sep=",")
# 
# head(canine_cts)
# dim(canine_cts)
