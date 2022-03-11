#Jenny Smith 

#February 16th, 2017

#purpose: Convert RNA-seq Expression files into GSEA broad institute format for CBF-GLIS patients GSEA. 
#GSEA requires that samples be normalized for sample-to-sample comparisions, so input should be TPM, RPKM, or batch corrected reads.

cls <- function(groupA,groupB, fileName){
  library(magrittr)
  #groupA and GroupB are character vectors with the patients IDs in each group
  #filename is character vector
  N <- length(groupA) + length(groupB)
  g1 <- as.character(substitute(groupA))
  g2 <- as.character(substitute(groupB))
  
  line1 <- paste(N, "2", "1", sep="\t")
  line2 <- paste("#", g1, g2, sep="\t")
  line3 <- c(rep(g1, length(groupA)), rep(g2, length(groupB)))
  
  cat(c(line1,line2), sep="\n", file=fileName)
  cat(line3, sep="\t", file=fileName, append=TRUE)
}



gct <- function(df, groupA, groupB, fileName){
  #groupA and GroupB are character vectors with the patients IDs in each group
  #df is the expression matrix. genes as rownames. MUST be filtered to remove genes with low counts
  samples <- length(groupA) + length(groupB)
  genes <- nrow(df)
  
  line1 <- "#1.2"
  line2 <- paste(genes, samples, sep = "\t")
  
  
  NAME <- rownames(df)
  Description <- rep("na", nrow(df))
  
  df <- df[,c(groupA, groupB)]
  df <- cbind(Description, df)
  df <- cbind(NAME, df)
  
  cat(c(line1, line2), sep="\n", file = fileName)
  write.table(df, sep="\t", file=fileName, quote=FALSE, row.names = FALSE, append=TRUE)
}




