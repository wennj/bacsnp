#' SNP specificity determination
#'
#' @param vcfdf A dataframe with SNP frequencies. Data was created with vcf_transformation.
#' @param isolates A character vector for which the SNP specificity determination is performed.
#' @param which.rel Specifies the column of relativ SNP frequencies for which the specificity determination should be performed on.
#' @return A list of two elements. First element is the input dataframe with SNP specificities added. The second element is a counting talbe with all SNP specificity combinations.
#' @examples
#' data(bac.f)
#' bacsnp.specificity(bac.f, c("iIso1", "iIso2"), which.rel = "REL.ALT1")
#' @importFrom stats kmeans
#' @export
bacsnp.specificity <- function(vcfdf, isolates, which.rel = "REL.ALT1"){

  extract.ISO_information <- function(df, isoname){

    #which isolate is needed?
    iso <- isoname

    #get only the desired isolate dataframe:
    iso.df <- df[df$ISO %in% iso,]

    #but I added previously the SNP positions from the paper?!
    #those are not truely detected SNPs. They should be removed again:
    #iso.df <- iso.df[complete.cases(iso.df[,"REF"]),]
    dim(iso.df)

    #give me positions only where the REL.ALT1 is UNEQUAL to ZERO (> 0 % Frequency)
    iso.df <- iso.df[iso.df$REL.ALT1 != 0,]

    return(iso.df)

  }

  get.kmeanCluster_ISO <- function(mydf, whichREL, ncluster){

    #mydf <- SNP.kmeans.comb
    #whichREL <- "REL.ALT1"
    #ncluster <- 1

    #which column do we consider for the kmeans?
    col2nd <- whichREL
    pos <- "POS"

    #create a vector which to names of the two columns
    twocols <- c(pos, col2nd)

    #extracht the two columns from the dataframe:
    iso.pos.ref <- mydf[,twocols]
    #head(iso.pos.ref)

    #isoCluster <- stats::kmeans(iso.pos.ref[,col2nd], ncluster, nstart = 10)
    isoCluster <- kmeans(iso.pos.ref[,col2nd], ncluster, nstart = 10)

    return(isoCluster)

  }

  add.groupingVariable <- function(mydf, group_data){

    #add the grouping variable to the data.frame
    mydf <- cbind(mydf, group_data)

    #rename the column of the grouping variable to "CLUS"
    n <- colnames(mydf)
    n[which(n == "group_data")] <- "CLUS"
    colnames(mydf) <- n

    return(mydf)


  }

  get.NTforALT <- function(ALT.df, position, alt123){

    p <- position
    pALT <- as.character(ALT.df[ALT.df$POS == p,]$ALT)
    pALT <- unlist(strsplit(pALT, "[,]"))[alt123]
    return(pALT)
  }

  #####
  #vcfdf <- bac.f
  #isolates <- unique(bac.f$ISO)
  #which.rel <- "REL.ALT1"
  ####

  ####run the de novo detection

  #remove non informative NA sites:
  vcfdf <- vcfdf[!is.na(vcfdf$REL.ALT1),]

  SNP.kmeans.comb <- extract.ISO_information(vcfdf, isolates)

  kcluster <-  get.kmeanCluster_ISO(mydf = SNP.kmeans.comb, whichREL = which.rel, ncluster = 1)
  SNP.kmeans.comb.clus <- add.groupingVariable(mydf = SNP.kmeans.comb, group_data = kcluster$cluster)
  kcluster_centers <- kcluster$centers
  c.grID <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "elven", "twelve")
  c.grID <- c(c.grID, as.character(13:100))

  ######################################################
  spec.pos <- list()
  spec.pos.rbind <- data.frame()

  c <- 1
  while (c != dim(kcluster_centers)[1]+1) {

    SNP.kmeans.comb.clus.c <- SNP.kmeans.comb.clus[SNP.kmeans.comb.clus$CLUS == c,]
    POS.M <- unique(SNP.kmeans.comb.clus.c$POS)
    SPEC <- c()
    GROUP.ID <- c()
    NT.SPEC <- c()

    p <- 1
    while(p != length(POS.M)+1){
      mc <- as.character(SNP.kmeans.comb.clus.c[SNP.kmeans.comb.clus.c$POS == POS.M[p],]$ISO)
      SPEC[p] <- paste(mc,collapse="_")
      lmc <- length(mc)
      GROUP.ID[p] <- c.grID[lmc]
      NT.SPEC[p] <- get.NTforALT(SNP.kmeans.comb.clus.c,POS.M[p],1)

      p <- p+1
    }

    spec.pos[[c]] <- as.data.frame(cbind(POS.M, SPEC, GROUP.ID, NT.SPEC))
    spec.pos[[c]]$SPEC <- gsub('i', '',spec.pos[[c]]$SPEC)

    x <- as.data.frame(cbind(POS.M, SPEC, GROUP.ID, NT.SPEC))
    spec.pos.rbind <- rbind(spec.pos.rbind,x)
    spec.pos.rbind$SPEC <- gsub('i', '',spec.pos.rbind$SPEC)

    c <- c+1

  }
  spec.pos.comb <- spec.pos.rbind
  ################################################END.


  combinations <- as.data.frame(table(spec.pos.comb$SPEC))
  #combinations <- combinations[order(-combinations["Freq"]),]
  combinations <- combinations[order(-combinations$Freq), ]

  spec.pos.all <- spec.pos.comb

  as.data.frame(table(spec.pos.all$SPEC))

  knownPOSdf <- spec.pos.all

  ######################################################
  add.knownSNPpos <- function(vcfdf.covff){

    isolates <- as.vector(unique(vcfdf.covff$ISO))

    df <- data.frame()
    i <- 1
    while(i!=length(isolates)+1){
      isolate <- vcfdf.covff[vcfdf.covff$ISO == isolates[i],]
      isolate <- merge(isolate, knownPOSdf, by.x = "POS", by.y = "POS.M", all = TRUE)
      isolates[i]

      ############change from > to >=

      #ask if there are rows in which the isolate position in the dataframe is NA:
      if(dim(isolate[is.na(isolate$ISO) == TRUE,])[1] >= 1){
        #if this is the case than ...
        #set the total COV to minCOV +1:
        #isolate[is.na(isolate$ISO) == TRUE,]$COV <- sum(minCOV, minALTcov,1)
        isolate[is.na(isolate$ISO) == TRUE,]$COV <- 100
        #set the coverage of reference (COV.REF) to minCOV + 1:
        isolate[is.na(isolate$ISO) == TRUE,]$COV.REF <- 100


        #set the coverage of all alternative coverages (COV.ALTx) to ZERO:

        #identify the number of alternatives found
        alternatives <- c("REL.ALT1","REL.ALT2","REL.ALT3")
        noALT <- length(which(colnames(vcfdf) %in% alternatives))
        #dpr <- geno(vcf)$DPR
        #noALT <- max(lengths(dpr))-1

        ##create a vector with all alternative coverage names:
        strCOV.ALT <- c()
        j <- 1
        while(j!=noALT+1){
          strCOV.ALT[j] <- paste("COV.ALT", j, sep = "")
          j <- j+1
        }

        #set the coverage of all alternative coverages (COV.ALTx) to ZERO:
        isolate[is.na(isolate$ISO) == TRUE,][strCOV.ALT] <- 0

        #put the isolate name in ISO:
        isolate[is.na(isolate$ISO) == TRUE,]$ISO <- isolates[i]

      }

      isolate$SPEC <- as.character(isolate$SPEC)

      if(length(isolate[is.na(isolate$SPEC) == TRUE,]$SPEC) != 0){
        isolate[is.na(isolate$SPEC) == TRUE,]$SPEC <- as.character(isolates[i])
      }

      isolate$GROUP.ID <- as.character(isolate$GROUP.ID)

      if(length(isolate[is.na(isolate$GROUP.ID) == TRUE,]$GROUP.ID)!=0){
        isolate[is.na(isolate$GROUP.ID) == TRUE,]$GROUP.ID <- "new"
      }

      df <- rbind(df, isolate)

      df$SPEC <- gsub('i', '',df$SPEC)

      i <- i+1
    }

    df$POS <- as.numeric(df$POS)
    df <- df[order(df$ISO, df$POS),]
    rownames(df) <- NULL

    return(df)

  }

  data <- add.knownSNPpos(vcfdf)

  vcfdf.spec <- list(data, combinations)
  names(vcfdf.spec) <- c("data","combinations")

  return(vcfdf.spec)
}
