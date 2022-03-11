#Jenny Smith

library(dplyr)
library(readr)
library(stringr)
library(magrittr)

CDEs <- readRDS(file.path(CDE,"1031/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_9.3.19.RDS"))

TPMs <- read_csv(file.path(HOME,
                           "0000.00.03_Expression_Matrices/TARGET_AML_0531_1031_Ribodepleted_RNAseq_dupGenesRemoved_TPM.csv"))
head(TPMs[,1:5])
dim(TPMs[,1:5])

cts <- read_csv(file.path(HOME,
                          "0000.00.03_Expression_Matrices/TARGET_AML_0531_1031_Ribodepleted_RNAseq_dupGenesRemoved_FractionalCounts.csv"))
head(cts[,1:5])
dim(cts)

cts.grch38 <- readRDS(file = file.path(HOME,
                                       "0000.00.03_Expression_Matrices/TARGET_AML_RBD_Kallisto_Quant_GeneLevel_scaledTPM_counts.RDS"))

head(cts.grch38[,1:5])
dim(cts.grch38)


TPMs.grch38 <- readRDS(file= file.path(HOME,
                                       "0000.00.03_Expression_Matrices/TARGET_AML_RBD_Kallisto_Quant_GeneLevel_TPM.RDS"))
head(TPMs.grch38[,1:5])
dim(TPMs.grch38)


##### Subset CDEs

CDEs <- CDEs %>% 
  filter(!is.na(USI), Study =="AAML1031") %>% 
  select(USI,Study,`Age in years`,`WBC (x10^3/MicroLiter) levels`,
         `Bone marrow leukemic blast percentage (%)`,`Peripheral blasts (%)`,
         `WHO Classification (final: path then study entry)`,
         `Primary Fusion/CNV`,`Additional Fusions/CNV`,
         `NPM mutation?`, `CEBPA mutation?`,`WT1 mutation?`,
         `FLT3/ITD positive?`,`FLT3/ITD allelic ratio`) %>% 
  set_rownames(.$USI)

head(CDEs[,1:5])
dim(CDEs)

##Subset Expn Data for 1031 only
subset_1031 <- function(CDE,ExpnData,NBM=FALSE){
  #https://stackoverflow.com/questions/30604107/r-conditional-evaluation-when-using-the-pipe-operator
  if(!NBM){
    rm.samples <- "Kas|MV4|MPN[0-9]|Sort|BM[0-9]|RO[0-9]"
    
  }else{
    rm.samples <- "Kas|MV4|MPN[0-9]|Sort"
  }
  print(NBM)
  
  rm.samples <- grep(rm.samples,colnames(ExpnData), ignore.case = TRUE)
  if(length(rm.samples)>0){ExpnData <- ExpnData[,-rm.samples]}
  if(!any(grepl("[A-Za-z]", rownames(ExpnData)))){colnames(ExpnData)[1] <- "Gene_Name"}
  
  if(any(grepl("\\.", colnames(ExpnData)))){
    USI <- str_split_fixed(colnames(ExpnData),pattern = "\\.", n=5)[,3]
  }else{
    USI <- colnames(ExpnData)
  }
  
  colNames <- c("Gene_Name",intersect(USI, CDEs$USI))%>% 
    {if(NBM) c(., grep("BM[0-9]|RO[0-9]",colnames(ExpnData), value=TRUE)) else .} %>%
    paste(.,  collapse = "|")
  
  CDE <- CDE[grep(colNames, CDE$USI, value = TRUE),]
  ExpnData <- ExpnData[,grep(colNames, colnames(ExpnData), value = T)]
  
  return(list(ExpnData, CDE))
}


cts.subset <-  subset_1031(CDE = CDEs, ExpnData = cts)
dim(cts.subset[[1]]) #51573  1077
dim(cts.subset[[2]]) #1063   14 (OK 1063+13 replicates + 1 Gene_Name column)
table(grepl("rep",colnames(cts.subset[[1]])))
# FALSE  TRUE 
# 1064    13 
head(cts.subset[[1]][,1:5])


cts.grch38.sub <- subset_1031(CDE=CDEs, ExpnData = cts.grch38,NBM=TRUE)
cts.grch38.sub[[1]] <- cts.grch38.sub[[1]][grep("^ENSG", rownames(cts.grch38.sub[[1]])),]
dim(cts.grch38.sub[[1]])
dim(cts.grch38.sub[[2]])
table(grepl("rep",colnames(cts.grch38.sub[[1]])))
# FALSE 
# 1130 
head(cts.grch38.sub[[1]][,1:5])


TPMs.grch38.sub <- subset_1031(CDE=CDEs, ExpnData = TPMs.grch38,NBM=TRUE)
TPMs.grch38.sub[[1]] <- TPMs.grch38.sub[[1]][grep("^ENSG", rownames(TPMs.grch38.sub[[1]])),]
dim(TPMs.grch38.sub[[1]])
dim(TPMs.grch38.sub[[2]])
table(grepl("rep",colnames(TPMs.grch38.sub[[1]])))
# FALSE 
# 1130 
head(TPMs.grch38.sub[[1]][,1:5])






### Save the Data 

# write.csv(TPMs.2, 
#           file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_DiagnosticSamples_TPM.csv"), 
#           row.names = FALSE)
# write.csv(cts.subset[[1]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_DiagnosticSamples_FractionalCounts.csv"),
#           row.names = FALSE)
# 

write.csv(cts.grch38.sub[[1]],
          file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_DiagnosticSamples_FractionalCounts.csv"),
          row.names = FALSE)

write.csv(cts.grch38.sub[[2]],
          file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_ClinicalData.csv"),
          row.names = FALSE)







