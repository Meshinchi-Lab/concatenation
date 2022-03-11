#Jenny Smith 
#April 5, 2019 
#Purpose: This will search for specific fusions evidence in FOI (RNAseq fusion data) and ISCN strings. 

# Limitation: the very complex regex for the ISCNs need to be visually checked. there are likely some rare combinations of each translocation that I was not aware of.
#in addition it WILL MISS any three/mutli-way translocations where the two chromosomes are not back to back - like t(11;19) is found vs t(11;14;19) is missed for KMT2As 



classify_fusions <- function(FOI, ISCN){
  require(dplyr)
  require(tibble)
  
  ref <- data.frame(matrix(ncol = 3, nrow = 0,
                           dimnames = list(NULL,c("Fusion","Regex","Cyto.Regex")))) %>% 
    add_row(Fusion=c("NUP98-NSD1",
                     "NUP98-KDM5A",
                     "CBFA2T3-GLIS2",
                     "KMT2A-MLLT3", 
                     "KMT2A-MLLT10",
                     "KMT2A-MLLT4",
                     "KMT2A-ELL",
                     "KMT2A-MLLT1", 
                     "RUNX1-RUNX1T1", 
                     "CBFB-MYH11"), 
            Regex= c("NUP98-NSD1|NSD1-NUP98", 
                     "NUP98-KDM5A|KDM5A-NUP98", 
                     "CBFA2T3-GLIS2|GLIS2-CBFA2T3", 
                     "KMT2A-MLLT3|MLLT3-KMT2A", 
                     "KMT2A-MLLT10|MLLT10-KMT2A", 
                     "KMT2A-MLLT4|MLLT4-KMT2A", 
                     "KMT2A-ELL|ELL-KMT2A", 
                     "KMT2A-MLLT1(;|$)|MLLT1-KMT2A(;|$)", 
                     "RUNX1-RUNX1T1|RUNX1T1-RUNX1", 
                     "CBFB-MYH11|MYH11-CBFB"), 
            Cyto.Regex=c("t\\(\\d{0,2};?5;11;?\\d{0,2}\\)\\([pq]?.{0,};?p35;p15.{0,}\\)", 
                         "t\\(\\d{02};?11;12;?\\d{0,2}\\)\\([pq]?.{0,};?p15;p13.{0,}\\)", 
                         "inv\\(16\\)\\(p13;?q24\\)", 
                         
                         "t\\(\\d{0,2};?9;11;?\\d{0,2}\\)\\s?\\([pq]?.{0,};?p2\\d;q23.{0,}\\)", 
                         "(t|dic|ins)\\(\\d{0,2};?10;11;?\\d{0,2}\\)\\([pq]?.{0,};?p1\\d;q23.{0,}\\)",  #46,XY,ins(10;11)(p12;q23q13)[18]/46,XY[2]
                         "t\\(\\d{0,2};?6;11;?\\d{0,2}\\)\\([pq]?.{0,};?q2\\d;q23.{0,}\\)", 
                         
                         "t\\(\\d{0,2};?11;19;?\\d{0,2}\\)\\([pq]?.{0,};?q23;?p13.1.{0,}\\)", 
                         "t\\(\\d{0,2};?11;19;?\\d{0,2}\\)\\([pq]?.{0,};?q23;?p13.[23].{0,}\\)", 
                         "t\\(\\d{0,2};?8;21;?\\d{0,2}\\)\\(q22;?q22\\)|ins\\(21;8\\)|ins\\(8;21\\)", 
                         "inv\\(16\\)\\(p13.{0,3}q22\\)|t\\(16;16\\)\\(p13.{0,3}q22|ins\\(16;16\\)\\(p13.1;?q22")) %>% 
    set_rownames(.$Fusion)
  
  res <- NULL
  for(r in 1:nrow(ref)){
    df <- ref[r, ]
    
    #Data Query
    FOI.Present <- any(grepl(df$Regex,FOI))
    Cyto.Present <- grepl(df$Cyto.Regex, ISCN)
    No.Cyto <- is.na(ISCN)
    
    #define classifications
    if(FOI.Present & No.Cyto){
      res <- c(res,"RNA seq only; no cyto avail")
    }else if(FOI.Present & Cyto.Present){
      res <- c(res,"both confirms")
    } else if(FOI.Present & !Cyto.Present){
      res <- c(res,"RNA seq only")
    } else if(!FOI.Present & Cyto.Present){
      res <- c(res,"cyto only")
    }else if(!FOI.Present & ! Cyto.Present){
      res <- c(res,NA)
    }
  }
  
  names(res) <- ref[,"Fusion"]
  if(sum(!is.na(res)) == 0){
    res <- NA 
    names(res) <- "none"
  }else if(sum(!is.na(res)) == 1 ){
    res <- res[! is.na(res)]
    
  }else if(sum(!is.na(res)) > 1){
    n <- paste(names(res[! is.na(res)]),collapse = "; ")
    res <- "fusion/cyto data conflict"
    names(res) <- n
  }
  
  res <- paste(names(res), res, sep=": ")
  return(res)
}