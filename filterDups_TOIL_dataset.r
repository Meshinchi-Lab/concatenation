#Jenny Smith 


#Jan 12, 2018 





#Purpose: filter duplicate genes from TOIL dataset with TCGA, TARGET, and Normal Tissues 

library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(reshape2)
library(tidyr)
library(tibble)


filterDups <- function(df, geneCol){
  #Function to input a expression matrix and remove duplicate gene entries
  #df is the expression matrix with a column for gene names as a character class. 
  #geneCol is the character vector with the column name. 
  #best usage: Arrange the dataframe by gene name BEFORE using this function. 
  
  library(genefilter)
  
  #Genes and duplicate genes in the dataframe
  genes <-unlist(df[,geneCol]) 
  dup.genes <- genes[duplicated(genes) | duplicated(genes, fromLast=TRUE)]
  
  #vector of all TRUEs to be the base vector to be updated based on rowVars. 
  keep <- rep(TRUE, nrow(df)) %>%
    set_names(genes)
  
  #subset the data frame for only duplicate genes
  df.dup <- df %>%
    filter_(paste(geneCol,"%in%","dup.genes"))
  
  #remove the larger dataframe from RAM
  rm(df)
  
  #For loop to determine which duplicate gene has the largest variance.
  for (dup in unique(dup.genes)){
    
    # print(dup)
    temp <- df.dup %>%
      filter_(paste(geneCol, "==", paste0('"', dup,'"'))) %>%
      # filter_(paste(geneCol, "==", paste0("\'",dup,"'"))) %>% #noticed possibly typo. fixed above
      select(-which(colnames(.) == geneCol))
    
    geneVars <- rowVars(temp)
    
    #if the variance is the same for all duplicate genes, pick the first one.
    if (max(geneVars) - min(geneVars) == 0){
      # print("equal var")
      update <- c(TRUE, rep(FALSE, nrow(temp)-1))
    }else{
      #else, select the gene with highest variation.
      update <- apply(temp, 1, function(x) var(x) >= max(geneVars))
      
      #more than one duplicate with the same max variance, select the first one.
      if (sum(update) > 1){
        update[which(update)]  <- c(TRUE,rep(FALSE, length(which(update))-1))
      }
    }
    
    keep[names(keep) == dup] <- update
  }
  
  
  return(keep)
}




# setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2018.01.02_BCL2_Expression/')
setwd('/fh/fast/meshinchi_s/workingDir/TARGET/NormalTissue_GTEX/RNA/mRNAseq/level3/gene/2016Sept_UCSC_Illumina_data/')

ID.Map <- read.table("/fh/fast/meshinchi_s/workingDir/TARGET/NormalTissue_GTEX/RNA/mRNAseq/metadata/gencode.v23.annotation.gene.probeMap", 
                     sep="\t", header = TRUE, stringsAsFactors = FALSE)

# head(ID.Map)

# *Unit is in log2 TPM + 0.001*
TCGA_Targ.log2 <- get(load("/fh/fast/meshinchi_s/workingDir/TARGET/NormalTissue_GTEX/RNA/mRNAseq/level3/gene/2016Sept_UCSC_Illumina_data/Rdata/TcgaTargetGtex_rsem_gene_tpm.RData"))
head(TCGA_Targ.log2[,1:5])
dim(TCGA_Targ.log2) #60498 by 19261 samples

options(scipen = 999)
TCGA_Targ <- TCGA_Targ.log2 %>%
  inner_join(., ID.Map[,1:2], by=c("sample"="id")) %>%
  select(gene, everything(), -sample) %>% #1,971 duplicate genes, no NAs
  mutate_if(is.numeric, function(x) 2^x)


options(scipen=999)
head(TCGA_Targ[,1:5])


TCGA_Targ.rmDups <- TCGA_Targ %>%
  arrange(gene) %>%
  mutate_if(is.numeric, function(x) ifelse(x > 0.001, x-0.001, round(x, digits = 3)-0.001)) %>%
  filter(filterDups(., geneCol = "gene"))


if (sum(duplicated(TCGA_Targ.rmDups$gene)) == 0){
  TCGA_Targ.rmDups <- TCGA_Targ.rmDups %>%
    column_to_rownames("gene")
}


head(TCGA_Targ.rmDups[,1:5])
dim(TCGA_Targ.rmDups[,1:5]) #58,531 genes 


save(TCGA_Targ.rmDups,file="Rdata/TcgaTargetGtex_rsem_dupGenesRemoved_tpm.RData" )
write.csv(TCGA_Targ.rmDups, file = "TcgaTargetGtex_rsem_dupGenesRemoved_tpm.csv")
# TCGA_Targ.rmDups <- get(load("/TcgaTargetGtex_rsem_geneSymbol_dupGenesRemoved_tpm.RData"))



# #Subset for those AML and whole blood 
# 
# toil <- read.table("/fh/fast/meshinchi_s/workingDir/TARGET/NormalTissue_GTEX/Clinical/TcgaTargetGTEX_phenotype.txt", 
#                    sep="\t", header = TRUE, stringsAsFactors=FALSE)
# 
# head(toil)
# # dim(toil) #19,131 by 7
# 
# sel <- toil %>%
#   filter(grepl("Blood|Myeloid Leukemia", primary.disease.or.tissue)) %>%
#   mutate(colname=gsub("\\-", ".",  sample)) %>%
#   select(sample, colname, everything())
# 
# 
# # sel
# dim(sel) #738 by 8
# table(sel$primary.disease.or.tissue, sel$X_study)
# 
# subset <- TCGA_Targ.rmDups %>%
#   select(intersect(sel$colname, colnames(TCGA_Targ.rmDups)))
# 
# 
# # dim(subset) #738 by 58,581
# head(subset[,1:5])
# 
# # save(subset,file= "TcgaTargetGtex_NormBlood_AML_rsem_geneSymbol_rmDups_tpm.RData")
# 
# 
# #The filter dups function works on smaller subsets of this dataset, but when run on the entire set,
# #pRNA is THEONLY one that still has 3 entries?? why
# #Becuase all 4 have variance == max variance. This has been updated now. 
# 
# TCGA_Targ.sub <- TCGA_Targ[,1:100] %>%
#   filter(grepl("^L|^M|^N|^O|^P", gene, ignore.case=TRUE)) %>%
#   mutate_if(is.numeric, function(x) ifelse(x > 0.001, x-0.001, round(x, digits = 3)-0.001)) %>%
#   filter(filterDups(., geneCol = "gene")) %>%
#   column_to_rownames("gene")
# 
# 
# 
# dim(TCGA_Targ.sub)
# head(TCGA_Targ.sub[,1:20])
# 
# 
# 
# t <- subset(TCGA_Targ.sub, TCGA_Targ.sub$gene == "pRNA")
# t
# 
# test1 <- filterDups(t, geneCol = "gene")
# sum(test1) #4 are true??? yes, because all 4 have variance == max variance. 
# which(test1)
# test1[which(test1)]  <- c(TRUE,rep(FALSE, length(which(test1))-1))
# test1
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
