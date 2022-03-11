#Jenny Smith

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(magrittr)
# setwd(file.path(SCRATCH,"jlsmith3"))
setwd(file.path(PROJHOME, "2021.02.09_lncRNA_Connerty"))



### Clinical Data

CDEs.old <- readRDS(file.path(CDE,"00_Archive/1031/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_9.3.19.RDS"))
CDEs.full <- readRDS(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_12.09.20.RDS")) 

inelig <- CDEs.full %>% 
  filter(USI != "Unknown") %>% 
  filter(Eligibility_Comments=="remove") %>% 
  pull(USI)

CDEs.full <- CDEs.full %>%  
  filter(USI != "Unknown", !is.na(USI)) 

dim(CDEs.full) #2314  150


### Manifest 

manifest <- read.csv(file.path(TARGET, "SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_02.04.21.csv"))

dim(manifest)
head(manifest) #3045   15

### Counts

# cts.grch38 <- readRDS(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/TARGET_AML_RBD_Dx_Rlps_NBM_MPN_Kallisto_Quant_GeneLevel_scaledTPM_counts.RDS"))
cts.grch38 <- readRDS(file = file.path(HOME,
                                       "0000.00.03_Expression_Matrices/TARGET_AML_RBD_Kallisto_Quant_GeneLevel_scaledTPM_counts.RDS"))
# write.csv(cts.grch38, "TARGET_AML_RBD_Kallisto_Quant_GeneLevel_scaledTPM_counts.csv")
# table(duplicated(colnames(cts.grch38))) #FALSE no duplicates
head(cts.grch38[,1:5])
dim(cts.grch38)

TPMs.grch38 <- readRDS(file.path(PROJHOME, "0000.00.03_ExpressionMatrices/Kallisto_GRCh38_Gencode_v29/TARGET_AML_RBD_Dx_Rlps_NBM_MPN_Kallisto_Quant_GeneLevel_Abundance_TPM.RDS"))
# TPMs.grch38 <- readRDS(file= file.path(HOME,
#                                        "0000.00.03_Expression_Matrices/TARGET_AML_RBD_Kallisto_Quant_GeneLevel_TPM.RDS"))
head(TPMs.grch38[,1:5])
dim(TPMs.grch38)


##### Subset CDEs

CDEs  <- CDEs.old %>% 
  inner_join(., select(CDEs.full,
                       Eligibility_Comments, USI,
                       Primary.Fusion, Primary.CNV,
                       Additional.Fusions.CNV,
                       WT1.mutation.,
                       ISCN_update=ISCN),
             by="USI") %>% 
  mutate(`Additional Fusions/CNV`=Additional.Fusions.CNV, 
         `WT1 mutation?`=WT1.mutation.,
         ISCN=ISCN_update) %>% 
  
  #Colnames changed overtime ... 
  filter(!is.na(USI), Study =="AAML1031") %>%
  select(USI,Study,
         Eligibility_Comments,
         `Age in years`,Sex,
         `WBC (x10^3/MicroLiter) levels`,
         `Bone marrow leukemic blast percentage (%)`,
         `Peripheral blasts (%)`,
         `WHO Classification (final: path then study entry)`,
         `Primary Fusion`=Primary.Fusion,
         `Primary CNV`=Primary.CNV,
         `Additional Fusions/CNV`,
         
         `NPM mutation?`, `CEBPA mutation?`,`WT1 mutation?`,
         `FLT3/ITD positive?`,`FLT3/ITD allelic ratio`,
          ISCN,
         `OS event ID`,`EFS event type ID`,
         `CR status at end of course 1`,`CR status at end of course 2`,
         `Year of Last Follow Up`, `OS time (days)`, `EFS time (days)`) %>%
  mutate(RNA_seq_TimePoint="diagnostic") %>%
  set_rownames(.$USI)

head(CDEs[,1:5])
dim(CDEs) #1175   26

# table(CDEs$Eligibility_Comments)
# table(CDEs$`EFS event type ID`)

CDEs.inelig <- CDEs %>% 
  filter(USI %in% inelig) %>% 
  select(USI, Study, Eligibility)

dim(CDEs.inelig) #84  3

  # filter(USI != "Unknown", Protocol=="AAML1031") %>% 
# write.csv(CDEs, "TARGET_AAML1031_ClinicalData.csv")



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
    print('Splitting Strings')
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


## AML counts 

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


## Normal Bone Marrow Counts 
NBM.samples <- manifest %>% 
  filter(grepl("NBM", Group)) %>% 
  filter(!grepl("replicate", Sample)) %>% 
  select(Sample, USI, Group, Lib_Prep)

# We already sent these to them....
# cts.grch38.NBM <- cts.grch38[grep("^ENSG", rownames(cts.grch38)),NBM.samples$Sample]
# 
# head(cts.grch38.NBM[,1:5])
# dim(cts.grch38.NBM)
# 
# TPMs.grch38.NBM <- TPMs.grch38[grep("^ENSG", rownames(cts.grch38)),NBM.samples$Sample]
# 
# head(TPMs.grch38.NBM[,1:5])
# dim(TPMs.grch38.NBM)



### Save the Data 

# write.csv(cts.grch38.NBM,"TARGET_Healthy_Controls_ribodepleted_Counts.csv")
# 
# write.csv(TPMs.grch38.NBM,"TARGET_Healthy_Controls_ribodepleted_TPM.csv")

# write.csv(CDEs,
#           "TARGET_AAML1031_ClinicalData_Update_2.22.21.csv",
#           row.names = FALSE)

# write.csv(NBM.samples, "TARGET_Healthy_Controls_ribodepleted_Samples.csv", row.names = F)

# write.csv(cts.grch38.sub[[1]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_ribodepleted_DiagnosticSamples_Counts.csv"),
#           row.names = TRUE)
# 
# write.csv(TPMs.grch38.sub[[1]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_ribodepleted_DiagnosticSamples_TPM.csv"),
#           row.names = TRUE)

# write.csv(cts.grch38.sub[[2]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AAML1031_ClinicalData.csv"),
#           row.names = FALSE)


## Results from Patrick 

files <- dir("Results/",full.names = T)

res <- lapply(files, read.csv)
names(res) <- dir("Results/")

#Duplicates bc ragged columns where "" is used as placeholder
# res$AC1204982.csv %>%  unlist() %>% duplicated() %>% table()
# FALSE  TRUE 
# 154    74
# lapply(res, head)   
# lapply(res, function(df) df %>%  unlist() %>% duplicated() %>% table()) 
              
res.long <- lapply(names(res), 
              function(x){
                    lncRNA <- gsub(".csv|-","", x)
                    col <- paste("Group",lncRNA, sep="_")
                    
                    res[[x]] %>%  
                      gather(!!col,USI) %>% 
                      filter(USI!="")
                  }) 
names(res.long) <- dir("Results/")

lapply(res.long, head)
lapply(res.long,dim)
#OK no duplicates
# dups <- lapply(res.long, function(x) filter(x, duplicated(USI)| duplicated(USI, fromLast = T)) %>% 
#                                               arrange(USI))


res_cat <- res.long %>%
  purrr::reduce(full_join, by = "USI") %>% 
  select(USI, everything())


head(res_cat)
dim(res_cat)
# table(res_cat$USI %in% inelig)
# write.csv(res_cat,"TARGET_AAAML1031_lncRNA_Groups_ChildrensCancerInstitute_and_Meshinchi.csv", row.names = F)



