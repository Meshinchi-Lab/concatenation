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
genome <-  "GRCh37"

if(genome=="GRCh37"){
  grch37_files <- dir(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/BCCA_GRCh37_Ensembl_v69"),
                      full.names=TRUE)
  # grch37_files
  
  ##Counts 
  grch37_cts_file <- grep("dupGenesRemoved_FractionalCounts", grch37_files,value=T)
  cts_grch37 <- readRDS(file.path(grch37_cts_file))
  
  gene_ids <- cts_grch37[,c(1:2)]
  cts_grch37 <- as.data.frame(cts_grch37)
  rownames(cts_grch37) <- cts_grch37$geneSymbol
  cts_grch37 <- cts_grch37[,-c(1:2)]
  
  
  ##TPM
  grch37_TPM_file <- grep("dupGenesRemoved_TPM", grch37_files, value = T)
  TPM_grch37 <- readRDS(file.path(grch37_TPM_file))
  
  gene_ids <- TPM_grch37[,c(1:2)]
  TPM_grch37 <- as.data.frame(TPM_grch37)
  rownames(TPM_grch37) <- TPM_grch37$geneSymbol
  TPM_grch37 <- TPM_grch37[,-c(1:2)]
  
  
  ## 0531 TPM 
  # polyA_files <-  dir(grch37_files[grep("PolyA", grch37_files)], full.names = TRUE)
  # TPM_0531_grch37 <- read.csv(file.path(grep("AAML0531_dupGenesRemoved_TPM", polyA_files, value=T)))
}


## Sample manifest and clinical data
sample_info <- read.csv(file.path(TARGET, "SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_08.12.21.csv")) 

dim(sample_info)


merged <- read.csv(file.path(CDE, "Merged/TARGET_AML_0531_1031_merged_CDEs_05.21.21.csv"), 
                   na.strings = c("N/A","#N/A","NA","^$", "^\\.$")) 

inelig <- merged %>% 
  filter(Eligibility_Comments == "remove") %>% 
  pull(USI)


## Subset the Sample Manifest
FUS.ERG_samples <- sample_info %>% 
  filter(Primary.Fusion=="FUS-ERG") 

dim(FUS.ERG_samples)
table(FUS.ERG_samples$Time_point)
table(FUS.ERG_samples$USI %in% inelig)

relapsed <- FUS.ERG_samples %>% 
  filter(Time_point=="relapse")

Regs <- merged %>% 
  filter(USI %in% relapsed$USI) %>% 
  select(Reg., Primary.Fusion, EFS.event.type.ID, SCT.in.1st.CR)

View(Regs)


#Save the Data 
# write.csv(Regs,file.path(SCRATCH,"jlsmith3/FUS-ERG_Relapse_Regs.csv"), row.names = FALSE)
# write.csv(FUS.ERG_samples,
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_GRCh37_FUS-ERG_RNA-seq_Manifest.csv"), row.names = FALSE)
# write.csv(cts_grch37[,FUS.ERG_samples$Sample],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_GRCh37_FUS-ERG_dupGenesRemoved_FractionalCounts.csv"), row.names = TRUE)

# write.csv(TPM_grch37[,FUS.ERG_samples$Sample],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_GRCh37_FUS-ERG_dupGenesRemoved_TPM.csv"), row.names = TRUE)


# head(TPM_grch37[,FUS.ERG_samples$Sample])
# head(cts_grch37[,FUS.ERG_samples$Sample])



