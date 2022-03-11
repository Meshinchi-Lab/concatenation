#Jenny Smith 

#June 26 2018 

#Purpose: Given an ensembl GTF file, parse it into a ID mapping file for transcript level data 



getIDmap <- function(GTF){
  #GTF is a dataframe from read.delim(gtf_file)
  library(dplyr)
  library(tibble)
  options(stringsAsFactors = FALSE)
  
  #standard ensembl GTF format and gencode GTF.
  tx <- GTF %>%
    filter(grepl("transcript", V3)) %>% 
    dplyr::pull(V9) %>% #use pull() to create vector from a single column 
    str_split(., pattern = "; ") %>% 
    lapply(., function(x) t(str_split(x, pattern = " ", simplify = TRUE))) %>% # can nest lappys to convert the entire entry into a dataframe with colname you want
    sapply(.,  function(x) set_colnames(x, value = x[1,])[-1,]) %>% #bapply ?
    sapply(., function(x) data.frame(as.list(x))) %>% 
    bind_rows(.)
  
  return(tx)
}

#Can use regex as well! 
# bioconductor: rtracklayer: import() will take a GTF in to GRanges object 
# plyrrange 
#ensemblDB / annotation hub - can download any version as GRanges or ensmblDB package 



# ##### LncRNAs from Gencode 
# 
# GTF <- read.delim("/fh/fast/meshinchi_s/workingDir/TARGET/Reference_Data/GRCh38/gtf/gencode.v29.long_noncoding_RNAs.gtf",
#                 comment.char="#", sep="\t", header = FALSE)
# 
# lncRNA <- getIDmap(GTF=GTF)
# 
# head(lncRNA)
# dim(lncRNA) # 29,566    17
# 
# write.csv(lncRNA,
#           "~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/gencode.v29.lncRNAs_Transcript.IDmap.csv",
#           row.names = FALSE)
# 
# table(lncRNA$transcript_type, useNA = "always")
# table(lncRNA$transcript_support_level,useNA = "always" )


###########

# GTF <- read.delim("/fh/fast/meshinchi_s/workingDir/TARGET/Reference_Data/GRCh38/gtf/gencode.v29.annotation.gtf",
#                 comment.char="#", sep="\t", header = FALSE)
# 
# all_genes <- getIDmap(GTF)
# write.csv(GTF, "~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/gencode.v29_Gene.ID.map.csv")


########### Ensembl v76

# GTF <- read.delim("/fh/fast/meshinchi_s/workingDir/TARGET/Reference_Data/GRCh38/Homo_sapiens.GRCh38.76.NoHeader.gtf",
# stringsAsFactors = FALSE, sep="\t", header=FALSE)

# head(GTF)
# v76 <- getIDmap(GTF=GTF)
# write.csv(v76, "~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/Homo_sapiens.GRCh38.76_Transcript.Gene.IDmap.csv")
# head(v76)
# dim(v76) #206,183  by    4
  
######################


############## ensembl v69



# GTF <- read.delim("/fh/fast/meshinchi_s/workingDir/TARGET/Reference_Data/GRCh37/gtf/Homo_sapiens.GRCh37.69.NoHeader.gtf",
# stringsAsFactors = FALSE, sep="\t", header=FALSE)
# head(GTF[,1:9])
# dim(GTF)
# 
# #No transcript of gene categories
# table(GTF$V3)
# 
# GTF[1:100,] %>%
#   filter(V3=="CDS")


# #standard ensembl GTF format and gencode GTF.
# tx <- GTF %>%
#     filter(grepl("CDS", V3)) %>%
#     dplyr::select(V9) %>%
#     unlist(.) %>%
#     str_split(., pattern = "; ") %>% 
#     lapply(., function(x) grep("gene_id|transcript_id|gene_name|transcript_name|protein_id", x, value=TRUE)) %>%
#     lapply(., function(x) t(data.frame(str_split(x, pattern = " ", simplify = TRUE)))) %>%
#     sapply(.,  function(x) set_colnames(x, value = x[1,])[-1,]) %>%
#     t(.) %>%
#     data.frame(.)
#   


# v69 <- getIDmap(GTF=GTF) #did this ever work?
# # write.csv(v76, "~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/Homo_sapiens.GRCh38.69_Transcript.Gene.IDmap.csv")
# head(v69)
# dim(v69) #206,183  by    4

#############################


# 
# GTF <- read.delim("/fh/fast/meshinchi_s/workingDir/TARGET/Reference_Data/GRCh37/gtf/Homo_sapiens.GRCh37.75.NoHeader.gtf", 
#                 stringsAsFactors = FALSE, sep="\t", header=FALSE)
# 
# 
# 
# head(GTF)
# 
# 
# # gene_id, transcript_id, gene_name, transcript_name
# 
# 
# tx <- GTF %>%
#   filter(grepl("transcript", V3)) %>%
#   select(V9) %>%
#   unlist(.) %>%
#   str_split(., pattern = "; ") %>% 
#   lapply(., function(x) grep("gene_id|transcript_id|gene_name|transcript_name", x, value=TRUE)) %>%
#   lapply(., function(x) t(data.frame(str_split(x, pattern = " ", simplify = TRUE),
#                                      stringsAsFactors=FALSE))) %>%
#   sapply(.,  function(x) set_colnames(x, value = x[1,])[-1,]) %>%
#   t(.) %>%
#   data.frame(.)
# 
# 
# head(tx)
# dim(tx) #215170      4
# 
# all(complete.cases(tx))

# tx[sample(nrow(tx), size = 5),]

# write.csv(tx, file = "~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/Homo_sapiens.GRCh37.75_Transcript.Gene.IDmap.csv", row.names = FALSE)


#Old Code from 5/13/19
# getIDmap <- function(GTF){
#   library(magrittr)
#   library(dplyr)
#   options(stringsAsFactors = FALSE)
#   
#   #standard ensembl GTF format and gencode GTF.
#   tx <- GTF %>%
#     filter(grepl("transcript", V3)) %>%
#     dplyr::select(V9) %>%
#     unlist(.) %>%
#     str_split(., pattern = "; ") %>% 
#     #added transcript support level - may be issue for GRCH37 which doesn't include this information
#     lapply(., function(x) grep("gene_id|transcript_id|gene_name|transcript_name|transcript_support_level", x, value=TRUE)) %>%
#     lapply(., function(x) t(data.frame(str_split(x, pattern = " ", simplify = TRUE)))) %>%
#     sapply(.,  function(x) set_colnames(x, value = x[1,])[-1,]) %>%
#     t(.) %>%
#     data.frame(.)
#   
#   return(tx)
# }



