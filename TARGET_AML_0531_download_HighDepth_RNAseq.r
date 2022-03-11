#Jenny Smith
#2/1/19
#purpose: download high-depth (discovery) RNA-seq data set for OHSU. 


#Set-up 

library(dplyr)
library(tibble)
library(rslurm)

options(stringsAsFactors = FALSE)

setwd("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/SequencingDataMatrix/")

#Read in the reference files

runTable <- read.delim("TARGET_AML_allSamples_SraRunTable.txt", sep="\t", header = TRUE)

head(runTable)
dim(runTable)


HD.USIs <- dir("/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/gene/2014Aug_BCCA_0531_Illumina_data/",
               pattern = "*.gene.*") %>%
  gsub(pattern = ".gene.quantification.txt", "", .) 


#Select Runs for downloading

HD.RNAseq <- runTable %>%
  filter(grepl("AML", study_name)) %>%
  filter(grepl("RNA", analyte_type)) %>%
  filter(Assay_Type == "RNA-Seq") %>%
  filter(Center_Name == "BCCAGSC") %>%
  filter(biospecimen_repository_sample_id %in% HD.USIs) %>%
  filter(AssemblyName=="") 
  # filter(duplicated(Sample_Name)| duplicated(Sample_Name, fromLast = TRUE))



head(HD.RNAseq)
dim(HD.RNAseq)

# write.csv(HD.RNAseq, "TARGET_AML_0531_Discovery_Cohort_RNAseq_SRA_RunTable.csv", row.names = FALSE)

#submit jobs the the cluster
head(HD.RNAseq$Run)
script <- "~/scripts/sbatch_jobs/sra_fastq-dump.sh"

jobs <- lapply(HD.RNAseq$Run, function(j) system(paste("sbatch", script, j),
                                                 intern = TRUE))

head(jobs)
length(jobs) #205 



#Identify Failed Jobs 

path="/fh/scratch/delete90/meshinchi_s/SRA/files/"

outfiles <- dir(path = path, pattern = "slurm.+out") %>%
  paste0(path, .)

#loop through the outfiles to find the errors 
failed=list()
i <- 1
for (file in outfiles){
  print(file)
  file <- readLines(file)
  idx <- grep("failed SRR[0-9]", file)
  
  if (length(idx) > 0){
    fail <- file[[idx]] %>%
      gsub("^.+(SRR[0-9].+)", "\\1", .)
    
    failed[[i]] <- fail
    i <- i + 1
  }
}

length(failed) #33 jobs

#restart failed jobs
jobs2 <- lapply(failed, function(j) system(paste("sbatch -e", paste0(j, ".stderr"),
                                                 "-o", paste0(j,".out"),
                                                 script, j),
                                                 intern = TRUE))










