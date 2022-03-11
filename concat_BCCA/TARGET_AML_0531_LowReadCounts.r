#Jenny Smith

#Feb 12, 2018
#Purpose: To examine the output of the multiQC report for Kallisto, which indicates low quality RNA-seq outliers. 



library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(stringr)

#setwd("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/metadata/2016Apr_BCCA_Illumina_QC/")
setwd(file.path(PROJHOME,"2017.01.05_RNAseq_BCCA28Apr2016_TPM_Conversion"))

samples <- read.csv("../2017July_BCCA_JSmith_Illumina_QC/TARGET_AML_1031_ReplicateSamples_forLastPlate_2.13.17.csv", stringsAsFactors = FALSE)
NBM <- read.csv("~/reference_mapping-files/Inventory_of_normal_bone_marrows.csv", stringsAsFactors = FALSE)
SeqSamples <- read.csv("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.08.18_RNAseq_TallyperPatient/TARGET_AML_AAML1031_NGS_Data_Availability.csv", stringsAsFactors = FALSE)



resTX <- read.table("resTx_data/mqc_kallisto_alignment_1.txt", sep="\t", header = TRUE)

dim(resTX)

quantile(resTX$Percent_Reads_Aligned)


presentIn1031 <- c("PASBHI", "PARAJX", "PASVYL", "PASVVS", "PARBIU")

resTX <- resTX %>%
  mutate(Total_Read_Count=Not.aligned + Pseudoaligned,
         Percent_Reads_Aligned=Pseudoaligned/Total_Read_Count *100,
         Mean_Read_Count=mean(Total_Read_Count)) %>%
  
  mutate(USI=str_split_fixed(Sample, pattern="-", n=5)[,3]) %>%
  filter(Percent_Reads_Aligned >= 65)


extra <- sample(resTX$USI, size=1)


extra %in% presentIn1031 #False


SeqSamples <- SeqSamples %>%
  mutate(USI=str_split_fixed(PATIENT_ID, "-", n=2)[,1])

NBMs <- setdiff(NBM$USI, SeqSamples$USI) %>%
  data.frame(V1=.,
             Plate=rep("NBMs", length(.)),
             stringsAsFactors = FALSE)



samples <- samples %>%
  bind_rows(c(Plate="Rep0531", V1=extra)) %>%
  bind_rows(NBMs)

getwd()
# write.csv(samples, "TARGET_AML_1031_ReplicateSamples_forLastPlate_2.13.17.csv", row.names = FALSE)


