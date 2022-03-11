#Jenny Smith

#9/1/2020

#Purpose: Check that AAML1031 RNA-seq is available on dbGaP. 

library(magrittr)
library(dplyr)
library(stringr)

dir(file.path(CDE,"Merged/"))
CDEs <- read.csv(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_7.08.20.csv"))

CDEs <- CDEs %>% 
  filter(!is.na(USI), USI != "Unknown")
dim(CDEs)



dir(file.path(TARGET,"SequencingDataMatrix/SRA/"))
runTable <- read.csv(file.path(TARGET,"SequencingDataMatrix/SRA/SraRunTable_TARGET_Cohorts_9.1.20.txt"))


# dim(runTable) #14785    73
head(runTable[,1:5])
View(head(runTable)) #biospecimen_repository_sample_id
# colnames(runTable)

table(runTable$Assay.Type)



TARGET_AML_RNAseq <- runTable %>% 
  filter(grepl("TARGET-20|TARGET-21",biospecimen_repository_sample_id), 
         study_name == "TARGET: Acute Myeloid Leukemia (AML)",
         Assay.Type == "RNA-Seq") %>% 
  mutate(USI=str_split_fixed(biospecimen_repository_sample_id, pattern = "-", n=5)[,3]) %>% 
  left_join(., select(CDEs, USI, Protocol), by="USI") %>%
  select(USI, biospecimen_repository_sample_id,Protocol, everything())




head(TARGET_AML_RNAseq[,1:10])
# tail(TARGET_AML_RNAseq[,1:10])
dim(TARGET_AML_RNAseq) #710

table(TARGET_AML_RNAseq$Protocol, useNA = 'ifany')

# write.csv(TARGET_AML_RNAseq,
#           file.path(TARGET, "SequencingDataMatrix/SRA/TARGET_AML_RNAseq_on_dbGaP_9.1.20.csv"),
#             row.names = FALSE)

