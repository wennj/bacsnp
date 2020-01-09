#' Transforming a VCF file into a dataframe.
#'
#' @param vcf A VCF file loaded with readVCF of the VariantAnnotation package.
#' @return A dataframe.
#' @import VariantAnnotation
#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors unstrsplit
#' @export
bacsnp.transformation <- function(vcf){

  beforeCOV.REF <- 7

  get.ALTs <- function(myvcf){

    dpr <- geno(myvcf)$DPR
    noALT <- max(lengths(dpr))-1

    #create a vector with all alternative coverage names:
    strPOS <- "POS"
    strCOV <- "COV"
    strCOV.REF <- "COV.REF"
    strREL.REF <- "REL.REF"
    strCOV.ALT <- c()
    strREL.ALT <- c()
    i <- 1
    while(i!=noALT+1){
      strCOV.ALT[i] <- paste("COV.ALT", i, sep = "")
      strREL.ALT[i] <- paste("REL.ALT", i, sep = "")
      i <- i+1
    }

    all.str <- c(strPOS, strCOV, strCOV.REF, strREL.REF, strCOV.ALT, strREL.ALT)

    return(all.str)
  }

  transform.VCFtolist <- function(myvcf, isolate){
    #prepare the parameters:

    #load the isolate names into a vector
    #isolate <- isolates
    #myvcf <- vcf

    #load the DPR data into an array
    dpr <- geno(myvcf)$DPR
    is.array(dpr)

    #load the DP data into an array
    dp <- geno(myvcf)$DP

    #load the DP data into an array
    gq <- geno(myvcf)$GQ

    #load the ALT data into an vector
    ar <- fixed(myvcf)$ALT

    #ar = IRanges::CharacterList(ar)
    #ar = S4Vectors::unstrsplit(ar, sep = ",")
    ar = CharacterList(ar)
    ar = unstrsplit(ar, sep = ",")


    #load the QUAL data into an vector
    qual <- fixed(myvcf)$QUAL

    #load the REF data into an vector
    rr <- fixed(myvcf)$REF
    rr = as.character(rr)

    #load the POS data into an vector
    p <- rownames(geno(myvcf)$GQ)
    p <- sub("[^:]*", "", p)
    p <- substring(p, 2)
    p <- strsplit(p, "_")
    p <- sapply(p, function(x) x[[1]])

    #what are the dimensions of our DPR-array (dpr)?
    dpr.dim <- dim(dpr)

    #What is the max amount of detected genotypes?
    #do it with lengthS not length!!
    gt.max <- max(lengths(dpr))

    #Some lists within the array are shorter. Let's fill them up to the maximum (gt.max).
    dprm <- lapply(dpr, `length<-`, gt.max)
    #replace all created "NA" with 0:
    dprm <- rapply(dprm, f=function(x) ifelse(is.na(x),0,x), how="replace" )
    #...and create a new array (matrix of lists)
    dprm <- matrix(dprm, nrow = dpr.dim[1], ncol = dpr.dim[2])
    colnames(dprm) <- isolate


    #prepare the result list argument:

    #POS	ISO	REF	ALT	COV	COV.REF	COV.ALT

    createALTnames <- function(nr){
      i <- 1
      altnames <- c("POS","QUAL","ISO","GQ","REF","ALT","COV","COV.REF",1:(nr-1))
      while(i!=nr){
        altnames[i+beforeCOV.REF+1] <- paste("COV.ALT", i, sep = "")
        i <- i+1
      }
      return(altnames)
    }

    getISOgt <- function(d, iso){

      d <- list()
      j <- 1

      while(j!=length(iso)+1){
        i <- 1
        c <- matrix(nrow = dpr.dim[1], ncol = length(createALTnames(gt.max)))
        colnames(c) <- createALTnames(gt.max)
        c[,"GQ"] <- unlist(matrix(lapply(gq[,j], function(x) x[[1]])))
        c[,"QUAL"] <- qual
        c[,"ISO"] <- paste("i",iso[j], sep = "")
        c[,"COV"] <- unlist(matrix(lapply(dp[,j], function(x) x[[1]])))
        c[,"ALT"] <- ar
        c[,"REF"] <- rr
        c[,"POS"] <- p

        while(i!=gt.max+1){
          c[,i+beforeCOV.REF] <- sapply(dprm[,j], function(x) x[[i]])
          i <- i+1
        }

        calc <- matrix(nrow = dim(c)[1], ncol = gt.max)
        k <- 1
        calc[,k] <- as.numeric(c[,beforeCOV.REF+k])/as.numeric(c[,"COV"])
        while(k!=gt.max){
          k <- k+1
          calc[,k] <- as.numeric(c[,beforeCOV.REF+k])/as.numeric(c[,"COV"])
        }

        n <- c(1:gt.max)
        l <- 1
        n[l] <- "REL.REF"
        while(l!=gt.max){
          l <- l+1
          n[l] <- paste("REL.ALT", l-1, sep = "")
        }
        colnames(calc) <- n

        c <- cbind(c,calc)

        d[[j]] <- c
        j <- j+1
      }
      names(d) <- iso

      m <- 1
      dm <- d[[m]]
      while(m!=length(d)){
        m <- m+1
        dm <- rbind(dm,d[[m]])
      }

      #setwd(dataDir.Input)

      return(d)
    }

    SNPisolate.list <- getISOgt(dprm,isolate)

    return(SNPisolate.list)

  }

  transform.VCFlist_to_VCFdf <- function(vcflist,alternatives,isolate){

    #how many isolates do we haven
    isolate <- paste("i", isolate, sep = "")


    all.str <- alternatives

    i <- 1
    vcfdf <- data.frame()
    while(i!=length(isolate)+1){
      vcfdf <- rbind(vcfdf, as.data.frame(vcflist[[i]]))
      i <- i+1
    }



    i <- 1
    while(i!=length(all.str)+1){
      vcfdf[,all.str[i]] <- as.numeric(levels(vcfdf[,all.str[i]])[vcfdf[,all.str[i]]])
      #as.numeric(vcfdf[,all.str[i]])
      i <- i+1
    }
    return(vcfdf)
  }

  get.all_ISO_names <- function(myvcf){
    s <- samples(header(myvcf))

    return(s)
  }

  isolates <- get.all_ISO_names(vcf)
  all.str <- get.ALTs(vcf)
  SNP.list <- transform.VCFtolist(vcf,isolates)
  SNP.df <- transform.VCFlist_to_VCFdf(SNP.list, all.str,isolates)

  SNP.df$QUAL <- as.numeric(as.character(SNP.df$QUAL))
  SNP.df$GQ <- as.numeric(as.character(SNP.df$GQ))

  #SNP.df[is.na(SNP.df)] <- 0

  SNP.df <- SNP.df[!is.na(SNP.df$REL.ALT1),]

  return(SNP.df)

}
