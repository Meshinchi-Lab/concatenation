library(dplyr)

# Define Functions 
source(file.path("~/scripts/conversion_scripts/Merge_Cat_FixDupIDs_Function.r"))

#Function for the TPM conversion. 
# Based on https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4
# Useage: sapply(cated$RPKM, RPKM_to_TPM)
RPKM_to_TPM <- function(RPKM){
  conversionFactor <- sum(RPKM) / 1E6
  TPM <- RPKM / conversionFactor
  
  return(TPM)
}


# Run 

out_path <- "/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.12.20_MYC-ITD_1031/Expression_Data"
path <- "/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/exon/2016Apr_BCCA_0531_Illumina_data/"


results <- dir(path, pattern="exon.normalized",full.names = TRUE)
length(results)

pattern="^.+(TARGET.+R)_.+" #Pattern to select the Target Barcode/GSC barcode

selected <- c(1,4,9,12) #Select the column indices
exon_data_0531 <- catExpnData(filenames = results, regex = pattern, cols = selected, header = FALSE)
names(exon_data_0531) <- c("gene_id","exon_number","counts","RPKM")

# str(exon_data_0531)
# lapply(exon_data_0531, function(x) head(x[,1:5]))

exon_data_0531[["TPM"]] <- apply(exon_data_0531$RPKM,2, RPKM_to_TPM)
table(apply(exon_data_0531[["TPM"]], 2, sum)) #all sum to 1 million

table(apply(exon_data_0531$gene_id,2,
            function(x) identical(x=x, y=exon_data_0531$gene_id[,1])))

exon_TPM_0531 <- exon_data_0531[["TPM"]] %>% 
  as.data.frame() %>% 
  mutate("gene_id"=exon_data_0531$gene_id[,1],
         "exon_number"=exon_data_0531$exon_number[,1]) %>% 
  select(matches("gene|exon"), everything())

exon_fractional_counts_0531 <- exon_data_0531[["counts"]] %>% 
  as.data.frame() %>% 
  mutate("gene_id"=exon_data_0531$gene_id[,1],
         "exon_number"=exon_data_0531$exon_number[,1]) %>% 
  select(matches("gene|exon"), everything())

exon_RPKM_0531 <- exon_data_0531[["RPKM"]] %>% 
  as.data.frame() %>% 
  mutate("gene_id"=exon_data_0531$gene_id[,1],
         "exon_number"=exon_data_0531$exon_number[,1]) %>% 
  select(matches("gene|exon"), everything())


# dim(exon_TPM_0531)


lapply(c("TPM","fractional_counts","RPKM"), function(x){
  write.csv(get(paste("exon",x,"0531",sep="_")), 
            file.path(out_path,paste0("TARGET_AML_low_depth_0531_exon_",x,".csv")))
})

