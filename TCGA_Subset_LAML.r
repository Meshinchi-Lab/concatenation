#Subset The LAML TCGA dataset


library(dplyr)


# setwd(dir = "/fh/fast/meshinchi_s/workingDir/TARGET/AML_TCGA/RNA/mRNAseq/level3/geneLevelExpn/BCCA_Illumina_data/")
setwd(file.path(WORKINGDIR, "TCGA/GDAC_LAML_TCGA/RNA/mRNAseq/level3/gene/2016.02.12_BCCA_Illumina_data"))


LAML <- read.csv("LAML_TCGA_TPM_geneExpression.csv", 
                 stringsAsFactors = FALSE, row.names = 1)


head(LAML[,1:5])
dim(LAML)



BCL2 <- LAML %>%
  filter(grepl("BCL2\\|",geneSymbol)) %>%
  mutate(Gene=str_split_fixed(geneSymbol, "\\|", n=2)[,1]) %>%
  select(geneSymbol, Gene, everything())



head(BCL2[,1:10])
dim(BCL2)


# write.csv(BCL2, "LAML_TCGA_BCL2_TPM_ExpressionMatrix.csv", row.names = FALSE)





