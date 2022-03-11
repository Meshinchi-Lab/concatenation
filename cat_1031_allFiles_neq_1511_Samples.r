#Jenny Smith


#Jenny Smith
#June 28, 2018 


#Purpose: to concatenate all the 1031 Kallisto RNAseq transcript abundances.  N=1,511 samples

library(dplyr)
library(magrittr)
library(stringr)
library(tidyr)
library(tibble)

source("~/scripts/conversion_scripts/Merge_Cat_FixDupIDs_Function.r")

addCols <- function(df,symbol,id){
  library(dplyr)
  
  df <- df %>%
    mutate(geneSymbol=symbol,
           gene_id=id) %>%
    select(geneSymbol, gene_id, everything())
  
  return(df)
}

setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.10.09_Concatenate_1031_RNAseq/')


filepath.1031 <-  "/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/transcript/2017July_BCCA_1031_Kallisto_Illumina_data/"

allfiles <- paste(filepath.1031,dir(path = filepath.1031, pattern = "^.+GRCh37.87.+.tsv$"),sep="")

target <- grep("Sorted|Kas|MV4", ignore.case=TRUE, allfiles,invert = TRUE, value=TRUE) #1,112 samples

length(target)

#Run Concatentation 

pattern <- "^.+\\/([BPRMK][A-Za-z0-9].+R)\\_.+"
selected.kallisto <- c("target_id",	"length",	"eff_length",	"est_counts",	"tpm")

cated <- catExpnData(filenames = target, regex = pattern, cols = selected.kallisto, header = TRUE)

lapply(cated, dim)  


cated.1031 <- lapply(names(cated)[2:5], function(x) addCols(df=as.data.frame(cated[[x]]), 
                                                                  symbol = rep(NA,nrow(cated[[x]])), 
                                                                  id=cated$target_id[,1])) %>%
  set_names(names(cated)[2:5])

lapply(names(cated.1031), function(x) write.csv(cated.1031[[x]][,-1],
                                                 paste0("TARGET_AML_AAML1031_Kallisto_Transcript_RNASeq_",x, ".csv"),row.names = FALSE))



