#Jenny Smith
#7/13/21

#Create a SOW manifest for BCCA for Canine AML

library(dplyr)
library(tibble)

setwd(file.path(HOME,"2021.07.13_CSU_Canine_AML_SOW"))

# Read in the Sample Info
canine_manifest <- file.path(PROJHOME,"2021.03.08_CSU_Canine_AML/2021 CSU canine acute leukemia samples.xlsx")
all_samples <- openxlsx::read.xlsx(canine_manifest, sheet="Masterlist")
canine_info <- openxlsx::read.xlsx(canine_manifest,sheet="Patient Info")
fqs <- dir(file.path(SCRATCH,"jlsmith3/CSU_Canine_AML/fastq"))
fqs <- data.frame(Sample.identifier=gsub("_.+","", grep("_1.fq", fqs, value=T)),
                  filenames=paste(grep("_1.fq", fqs, value=T), 
                                  grep("_2.fq",fqs, value=T), sep=","))

# Load the excel workbook
wb <- openxlsx::loadWorkbook("external_import_for_collab.xlsx")
names(wb)

# Read in the worksheets as data frames
Patient <- openxlsx::readWorkbook(wb, sheet="Patient", check.names = FALSE, sep.names = " ")

head(Patient)

Sample <- openxlsx::readWorkbook(wb, sheet="Sample", check.names = FALSE, sep.names = " ")

head(Sample)

Library <- openxlsx::readWorkbook(wb, sheet="Library", check.names = FALSE, sep.names = " ")

head(Library)

# Update the information 
Patient.update <- Patient %>% 
  add_row(Identifier=c(all_samples$CSU.Sample.number)) %>% 
  left_join(., select(canine_info, CSU.Sample.number, S=Sex),
            by=c("Identifier"="CSU.Sample.number")) %>% 
  mutate_at(vars(Sex), ~S) %>% 
  select(-S) %>% 
  rename_all(~gsub("\\."," ", .))

head(Patient.update)
dim(Patient.update)

Sample.update <- Sample %>% 

  add_row(Identifier=all_samples$CSU.Sample.number) %>% 
  mutate_at(vars(Patient.identifier), ~c(all_samples$CSU.Sample.number)) %>% 
  mutate_at(vars(Disease), ~gsub("Normal.+|Peripheral.+","",all_samples$Sample.type)) %>% 
  mutate_at(vars(Disease.status), ~case_when(
    Disease == "" ~ "Normal",
    TRUE ~ "Diseased")) %>%
  left_join(., select(canine_info,CSU.Sample.number,
                      Tentative.dx, Sample.type),
            by=c("Identifier"="CSU.Sample.number")) %>% 
  mutate_at(vars(Anatomic.site), ~case_when(
    !is.na(Sample.type) ~ Sample.type,
    is.na(Sample.type) ~ "Cell line")) %>% 
  mutate_at(vars(Type), ~case_when(
    grepl("Bone",Anatomic.site) ~ "bone marrow",
    grepl("blood",Anatomic.site) ~ "blood",
    TRUE ~ tolower(Anatomic.site))) %>%

  select(-Tentative.dx, -Sample.type) %>% 
  rename_all(~gsub("\\."," ", .))

table(Sample.update$Disease)
table(Sample.update$`Disease status`)
table(Sample.update$Type)


head(Sample.update)
dim(Sample.update)

Library.update <- Library %>% 
  add_row(Identifier=all_samples$CSU.Sample.number) %>% 
  mutate_at(vars(Sample.identifier), ~all_samples$RNAseq.identifier) %>% 
  mutate_at(vars(Nucleic.acid.type), ~"RNA") %>% 
  mutate_at(vars(Protocol), ~c("Poly A mRNA")) %>% 
  mutate_at(vars(`#.of.sequencing.runs.(required.if.FASTQ.submission)`), ~1) %>% 
  mutate_at(vars(`All.data.provided?`), ~ "Yes") %>% 
  left_join(., fqs, by=c("Sample.identifier")) %>% 
  mutate_at(vars(`FASTQ.file.name(s)`), ~filenames) %>% 
  select(-filenames) %>% 
  rename_all(~gsub("\\."," ", .))

head(Library.update)
dim(Library.update)

#Save the updated information
openxlsx::removeWorksheet(wb,
                          sheet = "Patient")
openxlsx::removeWorksheet(wb,
                          sheet = "Sample")
openxlsx::removeWorksheet(wb,
                          sheet = "Library")

openxlsx::addWorksheet(wb,
                       sheet = "Patient")
openxlsx::addWorksheet(wb,
                       sheet = "Sample")
openxlsx::addWorksheet(wb,
                       sheet = "Library")

openxlsx::writeData(wb,
                    sheet = "Patient",
                    Patient.update,
                    keepNA=FALSE,
                    rowNames=FALSE)
openxlsx::writeData(wb,
                    sheet = "Sample",
                    Sample.update,
                    keepNA=FALSE,
                    rowNames=FALSE)
openxlsx::writeData(wb,
                    sheet = "Library",
                    Library.update,
                    keepNA=FALSE,
                    rowNames=FALSE)
openxlsx::saveWorkbook(wb,
                       "Meshinchi_Canine_AML_external_import_for_collab_07.21.21.xlsx",
                       overwrite = T)
