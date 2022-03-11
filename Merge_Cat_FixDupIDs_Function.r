#Jenny Smith

#May 12 2017

#Purpose: To combine expression data and clinical data for analysis. 



fixDupIDs <- function(df, IDs=NULL,type){
  #use for dataframes with repetitive patient USIs which need to be set as rownames or column names. 
  #this happens with paired sample data. 
  #IDs is the column name containing the IDs which should be set as rownames. 
  #type is either "rownames", or "colnames"
  
  if (type == "rownames"){
    dups <- df[,IDs][which(duplicated(df[,IDs]))] 
    idx <- which(duplicated(df[,IDs]))
    for (i in 1:length(dups)){
      name <- paste(dups[i], ".1", sep="")
      df[idx[i], IDs] <- name
    }
    rownames(df) <- df[,IDs]
  }else if (type == "colnames"){
    cols <- colnames(df)
    idx <- which(duplicated(colnames(df)))
    for ( i in 1:length(idx)){
      name <- paste(cols[idx[i]], ".1", sep="")
      cols[idx[i]] <- name
    }
    colnames(df) <- cols
  }
  return(df)
}



#this needs to be TESTED
merge_CDE_Expn <- function(clinData, expnMatrix,geneList,unit="",ID.col="USI"){
  #expnData is a data frame with patient IDs as the column names, genes as rownames
  #Clindata is a data frame with patient IDs in rownames.  
  #geneList is a character vector  length 1 or greater
  #unit is a character vector length 1, which states the unit, eg RPKM, TPM, log2TPM
  #ID.col is the column in the clinData dataframe which will merge with the USIs from the expn matrix. 
  
  library(tibble)
  library(dplyr)
  
  #subset for genes of interest
  expnMatrix <- expnMatrix[rownames(expnMatrix) %in% geneList, ] 
  
  
  if(nrow(expn)==1){
    df <- data.frame(expn=unlist(expnMatrix)) %>%
      rownames_to_column() %>%
      set_colnames(c("Sample.ID",paste0(geneList,"_",unit)))
      
  }else{
    df <- expnMatrix %>%
      rownames_to_column("Gene") %>%
      gather(Sample.ID,val, -Gene) %>%
      spread(Gene,val) %>%
      rename_at(vars(2:NCOL(.)), funs(paste(., unit,sep="_")))
  }
  
  df <-   df %>%
    mutate(Group=case_when(
        grepl("Kas|MV4", Sample.ID) ~ "CellLine",
        grepl("BM[0-9]|RO[0-9]", Sample.ID) ~ "NBM",
        grepl("sorted", Sample.ID, ignore.case = TRUE) ~ "FlowSorted", 
        grepl("MPN[0-9]", Sample.ID) ~ "MPN", 
        TRUE ~ "AML")) %>%
    
    mutate(USI=gsub("TARGET.[02]0.", "", Sample.ID)) %>%
    mutate(USI=gsub(".0[39]A.+|.14A.+","", USI)) %>%
    mutate(USI=gsub(".Sorted.+|.Unsorted$|\\.[AD].+", "", USI)) %>% 
    
    left_join(., clinData, by=c("USI"=ID.col))
  

  return(df)
}


#Function for the TPM conversion. 
# Based on https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4
RPKM_to_TPM <- function(RPKM){
  #This is for use with one patient column at a time, so use apply() or sapply()
  conversionFactor <- sum(RPKM) / 1E6
  TPM <- RPKM / conversionFactor
  return(TPM)
}

#attempt to increase speed. but no not at all. 
catExpnData <- function(filenames,regex, cols, header=FALSE,skipLines=0){
  library(magrittr)
  #filenames is a character vector of all filenames. 
  #regex is a string with the pattern to extract the patient ID , eg "^.+(Kasumi|MV4)", from filenames 
  #cols is the character vector or numeric vector of the columns to select and concatenate. 
  
  extract_cols <-function(filename,cols,skipLines=0){
    
    
    if(all(skipLines > 0 & header)){
      print("Skipped")
      # aFile <- readLines(filename)[-1] #remove first line with extra info. 
      # aFile <- str_split_fixed(aFile, pattern = "\t",n = length(cols)) %>% #split into a matrix
      #   set_colnames(.[1,] ) %>%  #set colnames from the first line 
      #   .[-1, ] #remove the header row from matrix
      aFile <- read.delim(filename, sep="\t", header=header, as.is=TRUE,comment.char = "#",skip=skipLines)
    }else{
      aFile <- read.delim(filename, sep="\t", header=header, as.is=TRUE,comment.char = "#")
    }
    
    output <- list()
    for ( k in 1:length(cols)){
      colname <- cols[k]
      col <- aFile[,colname]
      output[[colname]] <- col
    }
    return(output)
  }
  
  combineColumns <- function(extract_cols.res,colname){
    sapply(extract_cols.res, '[[', colname)
  }
  
  
  IDs <- gsub(regex, "\\1", filenames)
  columns <- lapply(filenames,extract_cols,cols=cols, skipLines=skipLines) %>%
    set_names(IDs)
  
  catedMatrices <- lapply(cols, combineColumns, extract_cols.res=columns)  %>%
    set_names(cols)
  
  
  return(catedMatrices)
}


catRbind <- function(filenames,regex, header=FALSE, sep="\t",ID.Col.Name="Patient"){
  library(magrittr)
  library(data.table)
  #filenames is a character vector of all filenames. May need to paste() file paths.  
  #regex is a string with the pattern to extract the patient ID , eg "^.+(Kasumi|MV4)", from filenames 

  
  open_files <-function(filename,ID){
    aFile <- read.delim(filename, sep=sep, header=header, as.is=TRUE)
    
    if (nrow(aFile) < 1){
      aFile[1,] <- NA
    }
    
    aFile[,ID.Col.Name] <- rep(ID, nrow(aFile))
    return(aFile)
  }
  
  IDs <- gsub(regex, "\\1", filenames)
  
  files <- mapply(open_files,filenames,IDs, SIMPLIFY = FALSE) %>%
    set_names(IDs)
  
  catedFiles <- rbindlist(files,fill=TRUE)

  return(catedFiles)
}


mergeMolCols <- function(col1,col2){
  data <- paste(col1,col2, sep=" ") %>%
    gsub("^ | $", "", .) %>%
    sapply(., unique) 
  
  #change all empty strings
  data[data==""] <- "Not evaluated"
  
  return(data)
}


collapseDuplicates <- function(df,ID.column,duplicate){
  #df is the datframe with multiple patient entries to collapse
  #ID column is the column to match the USI as a character vector. 
  #duplicate is the USI of the dups as character vector. Must be unique'd. 
  
  #Usage collapseDuplicates is:
  # pathology.rmDups <- pathology[!(pathology$Patient.ID %in% dups), ]
  # collapseDups <- bind_rows(lapply(dups, function(x) collapseDuplicates(df=pathology, ID.column="Patient.ID", duplicate=x)))
  # pathology.rmDups <- rbind(pathology.rmDups, collapseDups)
  
  idx <- which(df[, ID.column] == duplicate)
  cde <- df[idx,]
  
  if (length(unique(cde)) == 1){
    cde <- unique(cde)
  }else{
    #Examine eac column seperately
    for (i in 1:ncol(cde)){
      #if all identical, just unique it
      if (length(unique(cde[,i])) == 1){
        cde[1,i] <- unique(cde[,i])
      }else{
        #otherwise, collap
        cde[1,i] <- paste(unique(cde[,i]), collapse = ";")
      }
    }
  }  
  
  #update the clinical annotations with only the merged cde.
  cde <- cde[1,]
  
  
  return(cde)
}



###########older version ############### 


# #attempt to increase speed. but no not at all. 
# catExpnData <- function(filenames,regex, cols, header=FALSE,removeFirstLine=FALSE){
#   library(magrittr)
#   #filenames is a character vector of all filenames. 
#   #regex is a string with the pattern to extract the patient ID , eg "^.+(Kasumi|MV4)", from filenames 
#   #cols is the character vector or numeric vector of the columns to select and concatenate. 
#   
#   extract_cols <-function(filename,cols){
#     output <- list()
#     
#     aFile <- read.delim(filename, sep="\t", header=header, as.is=TRUE)
#     
#     for ( k in 1:length(cols)){
#       colname <- cols[k]
#       col <- aFile[,colname]
#       output[[colname]] <- col
#     }
#     return(output)
#   }
#   
#   combineColumns <- function(extract_cols.res,colname){
#     sapply(extract_cols.res, '[[', colname)
#   }
#   
#   
#   
#   IDs <- gsub(regex, "\\1", filenames)
#   
#   columns <- lapply(filenames,extract_cols,cols=cols) %>%
#     set_names(IDs)
#   
#   catedMatrices <- lapply(cols, combineColumns, extract_cols.res=columns)  %>%
#     set_names(cols)
#   
#   
#   return(catedMatrices)
# }


# catExpnData <- function(fileNames,regex, cols, header=FALSE){
#   library(rowr)
#   #filenames is a character vector of all filenames. 
#   #regex is a string with the 
#   #cols is the character vector of the columns to select and concatenate. 
#   
#   i=1
#   output <- rep(list(list()), length(cols))
#   
#   IDs <- NULL
#   for (F in fileNames) {
#     aFile <- read.delim(F, sep="\t", header=header, as.is=TRUE)
#     ID <- gsub(regex, "\\1", F)
#     IDs <- c(IDs, ID)
#     # print(IDs)
#     
#     for ( k in 1:length(cols)){
#       col <- aFile[,cols[k]]
#       output[[k]][[i]] <- col
#     }
#     i <- i + 1
#     print(F)
#   }
#   
#   output <- lapply(output, function(x) do.call(cbind.fill, c(x, list(fill=NA)))) #NA option to avoid losing any information.
#   output <- lapply(output, setNames, nm=IDs)
#   # output <- lapply(output, do.call, what=cbind)
#   names(output) <- cols
#   
#   return(output)
# }




# merge_CDE_Expn <- function(clinData, expnMatrix,geneList, phenoVector=NULL){
#   #expnData is a data frame with patient IDs as the column names, genes as rownames
#   #Clindata is a data frame with patient IDs in rownames.  
#   
#   #subset for genes of interest
#   expnMatrix <- expnMatrix[rownames(expnMatrix) %in% geneList, ] 
#   
#   #ensure # of rows == # cols of expression data for later merging
#   expnMatrix <- expnMatrix[,intersect(colnames(expnMatrix), rownames(clinData))]
#   
#   if (is.null(phenoVector)){
#     tmp <- data.frame(t(expnMatrix))
#   }else{
#     phenoVector <- phenoVector[intersect(names(phenoVector), rownames(clinData))]
#     expnMatrix <- expnMatrix[, names(phenoVector)]
#     tmp <- data.frame(t(expnMatrix),
#                       Status=phenoVector)
#   }
#   
#   #merge the clinical and expn data
#   srt_clinData <- transform(merge(clinData, tmp, by.y=0, by.x=0),
#                             row.names=Row.names,
#                             Row.names=NULL)
#   
#   return(srt_clinData)
# }



