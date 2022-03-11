#Jenny Smith

#Oct 17, 2017
#Purpose: To examine the output of the multiQC report for Kallisto, which indicates low quality RNA-seq outliers. 



library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(stringr)

setwd("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/metadata/2017July_BCCA_JSmith_Illumina_QC/")
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/ggplot_Themes_Function.r")
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/Heatmaps_Function.r")
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/clusterAnalysis_Function.r")



CDE.1031 <- read.csv("~/reference_mapping-files/TARGET_AML_1031_CDE_cleaned_1.31.2018.csv",
                     stringsAsFactors = FALSE, row.names = 1)

dim(CDE.1031)



CDE.0531 <- read.csv("~/reference_mapping-files/TARGET_AML_current_asof_june30_2016_UPDATED_CLEAN_1.31.18.csv",
                     stringsAsFactors = FALSE, row.names = 1)

dim(CDE.0531)


getwd()

#batch info for miRNAseq and mRNAseq
batchInfo <- read.csv("~/reference_mapping-files/AAML1031_Reg_USI_conversion_for_miRNAseq_batchInfo.csv", 
                      stringsAsFactors = FALSE)


head(batchInfo)

#FastQC results
fastqc <- read.table("multiQC/fastqc_data/general_stats_table.tsv", header = TRUE,
                     stringsAsFactors = FALSE, sep="\t")

head(fastqc)

#Kallisto results 
resTX <- read.table("multiQC/resTx_data/mqc_kallisto_alignment_1.txt", sep="\t", header = TRUE)


resTX <- resTX %>%
  mutate(Total_Read_Count=Not.aligned + Pseudoaligned,
         Percent_Reads_Aligned=Pseudoaligned/Total_Read_Count *100,
         Mean_Read_Count=mean(Total_Read_Count),
         USI=str_split_fixed(.$Sample, pattern="-", n=2)[,1]) %>%
  inner_join(., batchInfo, by="USI")


head(resTX)
dim(resTX)


goodSamples <- resTX %>%
  filter(Percent_Reads_Aligned >= 40) %>%
  filter(USI %in% CDE.1031$USI)

#Plate 14 was all 0531/03P1
batchInfo %>%
  filter(Seq.Plate == "plate 14")


tab <- table(goodSamples$Seq.Plate)

samples <- list()
for (plate in names(tab)){
  temp <- goodSamples %>%
    filter(Seq.Plate == plate) %>%
    filter(Total_Read_Count >= mean(Total_Read_Count))
  # print(dim(temp))
  samp <- sample(temp$USI, size = 1 )
  samples[[plate]] <- samp
}

samples <- t(data.frame(samples)) %>%
  as.data.frame(.,stringsAsFactors=FALSE) %>%
  rownames_to_column("Plate") 

# write.csv(samples, "TARGET_AML_1031_ReplicateSamples_forLastPlate_2.13.17.csv", row.names = FALSE)

# sum(samples$V1 %in% c(CDE.1031$USI, CDE.0531$TARGET.USI.1))





#Read in the failed bams names (failed to be converted into fastq, were corrupted)

failed <- read.table("stderr/picard/failed_picard_12.20.17/failed_picard_samples.txt",
                     sep = "\t", stringsAsFactors = FALSE)


failed <- failed %>%
  mutate(USI=str_split_fixed(V1,"\\-",n=2)[,1]) %>%
  mutate(USI=gsub(" ","", USI)) %>%
  inner_join(., batchInfo, by=c("USI"="USI.1"))


head(failed)
tail(failed)

table(failed$Seq.Plate) #all are from plate 11. 


#Read in the counts
raw <- read.csv("~/RNA_seq_Analysis/2017.10.06_FLT3-ITD_DEGs/ExpressionData/TARGET_AML_AAML1031_dupGenesRemoved_FractionalCounts.csv",
                stringsAsFactors = FALSE, row.names = 1)
colnames(raw) <- str_split_fixed(colnames(raw), "\\.", n=2)[,1]





###########################

#Samples that failed seq duplication 
quantile(fastqc$X..Duplicate.Reads)
# 0%       25%       50%       75%      100% 
# 2.022836 11.785831 13.145077 15.285257 57.430301 

failedDups <- fastqc$X..Duplicate.Reads >= 50

# write.csv(fastqc[failedDups,], "TARGET_AML_1031_FailedSequenceDuplication.csv")
# write.csv(fastqc, "TARGET_AML_1031_SummaryStats_FastQC.csv")



##############################

#Summary Stats
sum(rownames(batchInfo) %in% colnames(raw))

table(resTX$Seq.Plate)
# plate 1 plate 10 plate 11 plate 12 plate 13 plate 14  plate 2  plate 3  plate 4  plate 5  plate 6  plate 7  plate 8  plate 9 
# 84       70       37       87       71        6       85       79       77       73       75       75       83       84 

range(resTX$Total_Read_Count)
#4,452,214 199,385,010

range(resTX$Percent_Reads_Aligned)
# 3.62451 73.60999

quantile(resTX$Percent_Reads_Aligned)
# 0%      25%      50%      75%     100% 
# 3.62451 33.63743 35.92046 39.54400 73.60999 

boxplot(resTX$Total_Read_Count)




#######################
#Select Cut-offs 
q1 <- quantile(resTX$Total_Read_Count, probs = c(.01,.99))
# 1%       99% 
# 10280673 142870425 
  
q5 <- quantile(resTX$Total_Read_Count, probs = c(.05,.95))
# 5%       95% 
# 70814970 107068273



q1.Percent_Aligned <- quantile(resTX$Percent_Reads_Aligned, probs = c(.01,0.99))
# 1%      99% 
# 26.34393 58.95753

q5.Percent_Aligned <- quantile(resTX$Percent_Reads_Aligned, probs = c(.05,0.95))
# 5%      95%
# 30.38806 46.57287 



#Subse the top and bottom 5% of Total_Read_Count reads
lowHighReads.5Per <- resTX %>%  
  filter(Total_Read_Count <= q5[1] | Total_Read_Count >= q5[2]) %>%
  add_row(USI="Mean Read Count", Total_Read_Count=86896615)



# write.csv(prop.table(table(lowHighReads.5Per$Seq.Plate)),"Top5Percent_Lowest_Highest_RawReadCounts_ByPlate_PropTable.csv", row.names = FALSE)



# tiff("Top5Percent_Lowest_Highest_RawReadCounts_ByPlate.tiff", height = 20, width = 8, units = "in", res = 600)
ggplot(lowHighReads.5Per, aes(x=USI, y=Total_Read_Count, fill=Seq.Plate)) +
  geom_bar(stat="identity") + 
  theme_numX +
  coord_flip()
# dev.off()


#Subset the top and bottom 1% of Total_Read_Count read counts
lowHighReads.1Per <- resTX %>%  
  filter(Total_Read_Count <= q1[1] | Total_Read_Count >= q1[2]) %>%
  add_row(USI="mean read count", Total_Read_Count=86896615)

# write.csv(lowHighReads.1Per, "TARGET_AML_1031_HighVarianceReads.csv")


# tiff("Top1Percent_Lowest_Highest_RawReadCounts_ByPlate.tiff", height = 16, width = 8, units = "in", res = 600)
ggplot(lowHighReads.1Per, aes(x=USI, y=Total_Read_Count, fill=Seq.Plate)) +
  geom_bar(stat="identity") + 
  theme_numX +
  coord_flip()
# dev.off()



#subset the top and bottom 1% of aligned
lowAligned.1perc <- resTX %>%
  filter(Percent_Reads_Aligned <= q1.Percent_Aligned[1] ) 

# tiff("Bottom1Percent_lowest_Reads_Aligned_bySeqPlate.tiff", height = 16, width = 8, units = "in", res = 600)
ggplot(lowAligned.1perc, aes(x=USI, y=Percent_Reads_Aligned, fill=Seq.Plate)) +
  geom_bar(stat="identity") + 
  theme_numX +
  coord_flip()
# dev.off()


#subset the top and bottom 5%
lowAligned.5perc <- resTX %>%
  filter(Percent_Reads_Aligned <= q5.Percent_Aligned[1] ) 


prop.table(table(lowAligned.5perc$Seq.Plate))



# tiff("Bottom5Percent_lowest_Reads_Aligned_bySeqPlate.tiff", height = 16, width = 8, units = "in", res = 600)
ggplot(lowAligned.5perc, aes(x=USI, y=Percent_Reads_Aligned, fill=Seq.Plate)) +
  geom_bar(stat="identity") + 
  theme_numX +
  coord_flip() 
# dev.off()





#PCA with Plate info
noPlateInfo <- which(!rownames(batchInfo) %in% colnames(raw))
RawCountsOnly <- colnames(raw)[which(! colnames(raw) %in%  rownames(batchInfo))]


batch <- batchInfo[intersect(rownames(batchInfo),colnames(raw)),] 
pheno <- batch$Seq.Plate %>% set_names(batch$USI.1)
cts <- raw[,intersect(rownames(batchInfo),colnames(raw))]

PCA.analysis <- PCA(expnData = cts, phenovector = pheno,round = TRUE)


# tiff("1031_rawCounts_PCA.tiff", height = 10, width = 10, res=600, units = "in")
PCA.analysis$pca_plot
# dev.off()





