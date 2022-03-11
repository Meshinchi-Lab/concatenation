#Jenny Smith
#Aug 6, 2020

library(dplyr)
library(readr)
library(stringr)
library(magrittr)


setwd(file.path(PROJHOME,"2017.05.25_Inv16_Russell_Rockne"))

#### Request

# Hello Rhonda and Soheil,
# 
# after a long time we are pleased to update you on this project.
# We recently published a mathematical modeling paper using a mouse model of inv16 AML, 
# 
# https://cancerres.aacrjournals.org/content/early/2020/07/02/0008-5472.CAN-20-0354
# 
# and we are now ready to apply the approach to human data.
# 
# We recently revisited the data you sent us previously (see attached and below for a summary). We believe we have promising results and would like to collaborate with you on an analysis of this data.
# 
# In particular, in order to verify our preliminary results, we need *raw count data or anonymized BAM files* for the samples you sent us previously (attached). Can you provide that data? Also, do you have any additional inv16 or healthy samples?
#   
#   We look forward to hearing from you and to collaborate on this project.
# 
# Thank you
# 
# Russ, Ya-Huei, Guido, and Sergio

#### Counts/TPMS
cts <- readRDS(file.path(PROJHOME,
                         "0000.00.03_ExpressionMatrices/TARGET_AML_DSAML_MPN_NBM_Ribodepleted_dupGenesRemoved_Fractionalcounts.RDS"))
head(cts[,1:5])
dim(cts) #48230  2345



##### Subset CDEs
# CDEs <- readRDS(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_7.08.20.RDS"))
CDEs <- readRDS(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_12.09.20.RDS")) 
dim(CDEs)

USIs <- str_split_fixed(colnames(cts),pattern = "\\.", n=5)[,3]

samples_1 <-  xlsx::read.xlsx("TARGET_AML_Inv_16_RNA_Seq_Samples-May2017.xlsx", sheetIndex = 1)
head(samples_1)

samples_2 <- read.csv("TARGET_AML_COH_Sample_Manifest.csv")
head(samples_2)
table(samples_2$AML_Subtype) #84 NBM samples 

all_samples <- unique(c(samples_1$Sample.ID, samples_2$USI))
length(all_samples) #248 samples

# table(samples_1$Sample.ID %in% CDEs$USI) #20 missing data
# table(samples_2$USI %in% CDEs$USI) #84 missing bc they are NBM samples 

# Arm B: Standard Chemotherapy (ADE 10+3+5) + Bortezomib
# Arm A: Standard Chemotherapy (ADE 10+3+5)



CDE.subset <- CDEs %>%  
  filter(USI != "Unknown") %>% 
  filter(USI %in% samples_1$Sample.ID | 
           USI %in% samples_2$USI) %>% 
  select(USI,Protocol,
         # Eligibility_Comments,
         Age.in.years:Blast.percent..by.flow.,
         ISCN,FLT3.ITD.positive.:CEBPA.mutation., 
         AAML1031_Treatment.Arm=Treatment.Arm,
         GO.Treatment_Gemtuzumab_ozogamicin=GO.Treatment, 
         Primary.Fusion:Additional.Fusions.CNV,
         Primary.CNV) %>% 
  mutate_at(vars(AAML1031_Treatment.Arm), ~case_when(
    grepl("Arm A", .) ~ "Arm A: Standard Chemotherapy (ADE 10+3+5)",
    grepl("Arm B", .) ~ "Arm B: Standard Chemotherapy (ADE 10+3+5)",
    TRUE ~ .)) %>% 
  mutate_at(vars(GO.Treatment_Gemtuzumab_ozogamicin), ~case_when(
    .=="Unknown" & Protocol == "AAML1031" ~ Protocol, 
    TRUE ~ .))


dim(CDE.subset) #178  21
# View(CDE.subset)
# write.csv(CDE.subset, "TARGET_AML_inv16_Clinical_Data.csv", row.names = F)

# table(CDEs$Eligibility_Comments, useNA = 'ifany') #OK


#### File Manifest 

# manifest <- read.csv(file.path(TARGET,"SequencingDataMatrix/TARGET_AML_Ribodepleted_Master_Manifest_8.5.20.csv"), row.names = 1)

manifest <- read.csv(file.path(TARGET,"SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_10.08.20.csv"),
                     row.names = 1)

# table(samples_2$USI %in% manifest$USI)

manifest_forCollab <- manifest %>% 
  filter(!is.na(USI), USI !="Unknown") %>% 
  filter((grepl("CBFB-MYH11", Primary.Fusion.CNV) & Group=="AML" & Time_point =="diagnostic") 
         | Group == "NBM")  %>%
  select(Sample,USI,Group,AML_Subtype,Tissue,Time_point)


head(manifest_forCollab)
dim(manifest_forCollab) #210   5
table(manifest_forCollab$Group)
table(manifest_forCollab$Time_point)
table(manifest_forCollab$Tissue)
table(manifest_forCollab$AML_Subtype)


##Subset Expn Data
cts_forCollab <- cts[,manifest_forCollab$Sample]
head(cts_forCollab[,1:5])
dim(cts_forCollab) #48230   210


### Save the Data 


write.csv(cts_forCollab,
          file.path(SCRATCH,"jlsmith3/TARGET_AML_Ribodepleted_RNAseq_COH_dupGenesRemoved_FractionalCounts.csv"),
          row.names = TRUE)

write.csv(manifest_forCollab,
          file.path(SCRATCH,"jlsmith3/TARGET_AML_COH_Sample_Manifest.csv"),
          row.names = FALSE)







