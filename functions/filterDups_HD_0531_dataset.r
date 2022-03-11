#Jenny Smith 


#Jan 17, 2018 


#Purpose: Filter duplicates from high depth 0531 discovery data set. 

library(stringr)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)


#Fixing differences noted in the new code.... 
# source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/rmDupGenes_Function.r")

setwd("~/RNA_seq_Analysis/0000.00.03_Expression_Matrices/")



IDmap <- read.csv("~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/GeneSymbol_EnsemblID_Conversion_GRCh37.69_FromBCCA.csv", stringsAsFactors = FALSE)
head(IDmap)



tpm.0531.HD <- read.csv("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/gene/concat_matrices/2014Aug_BCCA_0531_Concatenated_Illumina_data/TARGET_AML_TPM_Aug2014.csv", 
                        stringsAsFactors = FALSE, row.names = 1)
rpkm.0531.HD <- read.csv("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/gene/concat_matrices/2014Aug_BCCA_0531_Concatenated_Illumina_data/TARGET_AML_RPKM_Aug2014.csv", 
                        stringsAsFactors = FALSE, row.names = 1)
cts.0531.HD <- read.csv("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/gene/concat_matrices/2014Aug_BCCA_0531_Concatenated_Illumina_data/TARGET_AML_rawCounts_Aug2014.csv", 
                        stringsAsFactors = FALSE, row.names = 1)




originals <- list(TPM=tpm.0531.HD, RPKM=rpkm.0531.HD, FractionalCounts=cts.0531.HD)
lapply(originals,function(x) head(x[,1:5]))
lapply(originals, dim)

#Add Gene names 
addGenes <- function(df, IDmap){
  df <-  df %>%
    inner_join( ., IDmap, by=c("ensemblSymbol"="gene_id")) %>%
    select(geneSymbol, everything()) %>%
    arrange(geneSymbol)
}


df.gene.names <- lapply(originals, addGenes, IDmap=IDmap)
lapply(df.gene.names, function(x) head(x[,1:5]))
lapply(df.gene.names, function(x) sum(is.na(x[,"geneSymbol"]))) #no missing entries in geneSymbol
lapply(df.gene.names, dim)

#write to files
# lapply(names(df.gene.names), function(x) 
#   write.csv(df.gene.names[[x]], paste0("TARGET_AML_withGeneNames", x, "_Aug2014.csv")))


############## Filter Duplicates #############


rmDupGenes <-  function(expnData, geneCol){
  library(genefilter)
  library(magrittr)
  #expnData is the matrix with a column for gene names.
  #geneCol a character vector with the name of the column for gene names.
  dups <- unique(expnData[which(duplicated(expnData[, geneCol])), geneCol])  ##MUST ADDRESS NAs!!!
  
  #add regex anchors to duplicate gene names
  # dup.regex <- unique(dups) %>% paste0("^", ., "$")#add regex anchors
  dup.regex <- paste0("^", dups, "$") 
  
  genes <- expnData[,geneCol]
  vars <- genefilter::rowVars(expnData[,-(which(colnames(expnData) == geneCol))])
  names(vars) <- genes
  
  tokeep <- function(dup,vars){
    geneVar <-  vars[grepl(dup, names(vars))]
    dup <- gsub("\\^|\\$", "", dup) #strip the anchors
    keep <- which(names(vars) == dup & vars == max(geneVar))
    #if there are ties with same max variation.
    if (length(keep) > 1){
      keep <- keep[1]
    }
    return(keep)
  }
  
  keep <- sapply(dup.regex,tokeep,vars=vars)
  sel.Dups <- expnData[keep, ]
  noDups <- expnData[which(! expnData[,geneCol] %in% dups), ]
  # print(sapply(list(sel.dup, noDups), dim))
  
  remDups <- rbind(noDups, sel.Dups) %>% set_rownames(., .[,geneCol])
  remDups <- remDups[,-(colnames(remDups) == geneCol)]
  
  list <- list(dups,dup.regex, vars, keep, remDups)
  names(list) <- c("dups", "dup.regex", "vars", "keep", "remDups")
  return(list)
}

df.rmDups <- lapply(df.gene.names, function(x) rmDupGenes(expnData=x[,-2], geneCol="geneSymbol"))
dim(df.rmDups$TPM$remDups)


lapply(df.rmDups, function(x) dim(x$remDups))
lapply(df.rmDups, function(x) head(x$remDups[,1:5]))

# lapply(names(df.rmDups), function(x) write.csv(df.rmDups[[x]]$remDups, paste0("TARGET_AML_dupGenesRemoved_", x, "_Aug2014.csv")))


######## Subset the Dataframes


#Subset the datasets for Dx and Relapse samples
filterSamp <- function(df, pattern){
  df <- df %>%
    select(grep(pattern, colnames(.)))
  return(df)
}

df.rmDups <- lapply(df.rmDups, function(x) x$remDups)
dx <- lapply(df.rmDups,filterSamp, pattern="09A|03A")

lapply(dx, function(x) head(x[,1:5]))
lapply(dx, dim) #160 samples
# lapply(names(dx), function(x) write.csv(dx[[x]], paste0("TARGET_AML_DxSamples_dupGenesRemoved_", x, "_Aug2014.csv")))


relapse <- lapply(df.rmDups, filterSamp, pattern="04A|40A")
lapply(relapse, function(x) head(x[,1:5]))
lapply(relapse, dim) #47 samples
# lapply(names(relapse), function(x) write.csv(relapse[[x]], paste0("TARGET_AML_RelapseSamples_dupGenesRemoved_", x, "_Aug2014.csv")))



##### Old code 



# #Filter the TPMs for duplicates
# filterManyDf <- function(df, IDmap){
#  df <-  df %>%
#     inner_join( ., IDmap, by=c("ensemblSymbol"="gene_id")) %>%
#     select(geneSymbol, everything(), -ensemblSymbol) %>%
#     arrange(geneSymbol) %>%
#     filter(filterDups(., geneCol="geneSymbol")) %>%
#     column_to_rownames("geneSymbol")
# }
# 



# rmDups <- lapply(originals, filterManyDf, IDmap=IDmap)
# names(rmDups) <- c("TPM", "RPKM", "cts")
# 
# lapply(rmDups, function(x) head(x[,1:5]))
# lapply(rmDups, dim)
# # lapply(names(rmDups), function(x) write.csv(rmDups[[x]], paste0("TARGET_AML_dupGenesRemoved_", x, "_Aug2014.csv")))







