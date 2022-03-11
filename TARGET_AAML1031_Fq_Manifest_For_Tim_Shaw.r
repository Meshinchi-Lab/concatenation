# Jenny Smith
# 7/13/21 

#Fastq Manifest for Tim Shaw for samples used in the interactive tSNE plot

library(dplyr)
library(stringr)

setwd(file.path(PROJHOME,"2021.07.13_Interactive_tSNE"))

# Input samples used by Tim Shaw
files <- read.delim("sample_name_rna_library_type.txt")

head(files)
dim(files) #2887    1

length(unique(files$SampleName))


# Fastq Manifest (out-of-date)
fq_manifest <- read.csv(file.path(TARGET,"SequencingDataMatrix/Fastq_manifests/TARGET_AML_RBD_PolyA_AWS_S3_Fastq_Rename_Log_11.18.20.csv")) %>% 
  rowwise() %>% 
  mutate(SampleName=gsub("^.+picard_fq2/\\(.+)$", "\\1", final_fastq_filename) %>%
           str_split_fixed(., pattern = "_", n=5) %>% 
           .[,1:2] %>% 
           paste(., collapse = "_")) %>% 
  ungroup() %>% 
  select(Sample:Lib_Prep,SampleName) %>% 
  distinct()



head(fq_manifest)
View(head(fq_manifest))
dim(fq_manifest) #2834    4

# fq_manifest %>% 
#   filter(grepl("CBFA2T3_GLIS2", Sample))

#RBD Counts Matrix manifest

sample_info <- read.csv(file.path(TARGET,
                                  "SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_06.09.21.csv")) %>%
  bind_rows(., read.csv(file.path(TARGET,
                                  "SequencingDataMatrix/BAM_manifests/TARGET_AML_polyA_RNAseq_Bam_Manifest_10.02.20.csv"))) %>% 
  filter(Protocol != "HD_NBM")

dim(sample_info)
head(sample_info)
table(sample_info$Batch)



# Merge in the fastq file manifest 
files_lib_prep <- files %>% 
  left_join(., fq_manifest, by=c("SampleName")) %>% 
  
  #The cord blood samples and some cell lines are missing due to being added later 
  #So fill in the NA values and I can update the full  manifests
  mutate_at(vars(Sample), ~case_when(
    is.na(.) ~ gsub("-","\\.", str_split_fixed(SampleName, pattern="_", n=2)[,1]),
    TRUE ~ .)) %>% 
  mutate_at(vars(Lib_Prep), ~case_when(
    is.na(.) ~ "RBD",
    !is.na(.) ~ gsub("RBS","RBD", .))) %>% 
  mutate_at(vars(Batch), ~case_when(
    is.na(.) ~ "rem2",
    TRUE ~ .))
  # left_join(., select(sample_info, Sample,USI,Protocol:Lib_Prep),
  #           by="Sample")


dim(files_lib_prep)
head(files_lib_prep)

any(duplicated(files_lib_prep$SampleName))
any(duplicated(files_lib_prep$Sample))
table(files_lib_prep$Lib_Prep)


# write.csv(files_lib_prep,"TARGET_AML_Fastq_Files_Library_Prep_Info.csv", row.names = FALSE)


