#Jenny Smith 
# 7/13/2018
#Purpose: "Concatenate 1031 Exon Coverage"




setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.10.09_Concatenate_1031_RNAseq/')



library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(reshape2)
getwd()



source("~/scripts/conversion_scripts/Merge_Cat_FixDupIDs_Function.r")



addCols <- function(df,symbol,exon){
  library(dplyr)
  
  df <- as.data.frame(df)
  
  df <- df %>%
    mutate(geneSymbol=symbol,
           exon=exon) %>%
    select(geneSymbol, exon, everything())
  
  return(df)
}





#Identify Files to Be Concatenated
filepath <-  "/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/exon/2017July_BCCA_1031_RSubread_JSmith_Illumina_data/"

allfiles <- dir(path = filepath,
                pattern = ".exon.cts")

head(allfiles)
length(allfiles) #1,117 files (no stella)



#target matrix will have NBM, AML, and untreated cell line samples ("D1" == day 1)
target <- paste0(filepath, grep("^[RBPS][A-Z0-9\\-]", allfiles, value=TRUE)) #1,111 samples

length(target) #1111




#Columns Description

# For exon.cts from Rsubread
# 
# Geneid	
# Chr	
# Start	
# End	
# Strand	
# Length	
# RO02505-09A-01R_withJunctionsOnGenome_dupsFlagged.bam

#note: this may cause issues due to all the last columns having a different colname.
#Plus these have an addition line at the top with commands. So, I will 



#Begin Concatenation 


#Pattern to select the Target Barcode
pattern <- "^.+\\/([BPR][A-Z0-9].+R)\\_.+"

#Select the column indices 
selected <- c(1:7)


cated <- catExpnData(filenames = target[1:10], regex = pattern, cols=selected, header = TRUE, removeFirstLine = TRUE)


sapply(cated, dim) 


gc()
# save(cated, file="TARGET_AML_1031_ExonLevel_Subread_RNAseq_Cated.RData")
# load("TARGET_AML_1031_ExonLevel_RNAseq_Cated.RData")



#Check that the Gene Ids are in the Same Order

# apply(cated$TARGET$`1`, MARGIN=2,FUN=identical, y=cated$TARGET$`1`[,1])
sapply(cated$TARGET, function(mat) all(apply(mat, MARGIN=2,FUN=identical, y=mat[,1])))







#Add columns for the Gene Symbol and Ensembl Symbol 
# withCols <- lapply(cated$TARGET[3:5], addCols, 
#                    symbol = cated$TARGET$geneSymbol[,1], exon = cated$TARGET$Exon[,1])
# 
# 
# lapply(withCols, function(x) head(x[,1:5]))
# 
# 
# 
# withCols.cells <- lapply(cated$Cells[3:5], 
#                          addCols, symbol = cated$Cells$geneSymbol[,1], exon = cated$TARGET$Exon[,1])
# 
# lapply(withCols.cells, function(x) head(x[,1:5]))



#Save the output






#Session Info
sessionInfo()
