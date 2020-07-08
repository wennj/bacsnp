#' Filtering SNP frequency dataframe
#'
#' @param vcfdf A dataframe with SNP frequencies. Data was created with vcf_transformation.
#' @param min.abs.cov The minimal absolute coverage (COV) that is required to consider SNP positions. All positions with a miniam absolute coverage below this value will be removed.
#' @param min.abs.alt An absolute value for the minimal number of occurances of an alternative nucleotide.
#' @param min.rel.alt An relative value for the minimal frequency of an alternative nucleotide.
#' @return A vcf dataframe that was filtered by a required minimal total coverage (COV), minimal coverage of alternative nucleotides (COV.ALT) or a minimal relative frequency of the first alterntative (REL.ALT1).
#' @examples
#' bacsnp.filter(bac, min.abs.cov = 100)
#' bacsnp.filter(bac, min.abs.cov = 100, min.abs.alt = 10)
#' bacsnp.filter(bac, min.abs.cov = 100, min.abs.alt = 10, min.rel.alt = 0.05)
#' bacsnp.filter(bac, min.rel.alt = 0.05)
#' @export
bacsnp.filter <- function(vcfdf, min.abs.cov = 0, min.abs.alt = 0, min.rel.alt = 0){

  #vcfdf <- bac.df
  #min.abs.cov = 100
  #min.abs.alt = 10
  #min.rel.alt = 0.05


  #get the number of alternatives
  #noALT <- count(grepl("REL.ALT", colnames(vcfdf)))
  noALT <- sum(grepl("REL.ALT", colnames(vcfdf)), na.rm = TRUE)

  #create a vector with all alternative coverage names:
  strCOV.REF <- "COV.REF"
  strCOV.ALT <- c()
  i <- 1
  while(i!=noALT+1){
    strCOV.ALT[i] <- paste("COV.ALT", i, sep = "")
    i <- i+1
  }
  COV.str <- c(strCOV.REF, strCOV.ALT)

  #create a vector with all alternative coverage names:
  strREL.REF <- "REL.REF"
  strREL.ALT <- c()
  i <- 1
  while(i!=noALT+1){
    strREL.ALT[i] <- paste("REL.ALT", i, sep = "")
    i <- i+1
  }
  REL.str <- c(strREL.REF, strREL.ALT)

  #Filter by absolute total coverage (COV)
  if(min.abs.cov > 0){
    vcfdf <- vcfdf[vcfdf$COV >= min.abs.cov,]
  }

  #Filter by too low alternative frequency
  if(min.abs.alt > 0){
    vcfdf[,COV.str] <- lapply(vcfdf[,COV.str], function(x) ifelse(x < min.abs.alt, 0,x))
  }

  #moved down to its own section:
  #Filter by too low alternative frequency
  #if(min.rel.alt > 0){
  #  vcfdf[,REL.str] <- lapply(vcfdf[,REL.str], function(x) ifelse(x < min.rel.alt, 0, x))
  #}

  #Recalculate the total COV:
  vcfdf[,"COV"] <- as.vector(rowSums(vcfdf[,COV.str]))

  #Recalculate the relative frequencies
  #first search for positions where COV is ZERO
  COV.isZERO.df <- vcfdf[vcfdf$COV == 0,]
  vcfdf <- vcfdf[vcfdf$COV != 0,]

  #divide by the COV
  vcfdf[,REL.str] <- vcfdf[,COV.str]/vcfdf$COV

  vcfdf <- rbind(vcfdf,COV.isZERO.df)

  #------------------------------Filtering by REL.ALT

  #Filter by too low alternative frequency
  if(min.rel.alt > 0){
    vcfdf[,REL.str] <- lapply(vcfdf[,REL.str], function(x) ifelse(x < min.rel.alt, 0, x))

    #vcfdf[,strREL.REF] <- 1-rowSums(vcfdf[,REL.str[which(REL.str != strREL.REF)]])
    vcfdf[,strREL.REF] <- 1-rowSums(vcfdf[,strREL.ALT])
    #recalculate the COV values
    vcfdf[,COV.str] <- vcfdf$COV * vcfdf[,REL.str]

  }


  return(vcfdf)

}
