#' Transforming a VCF file into a dataframe.
#'
#' @param vcf A VCF file loaded with read.vcfR of the vcfR package.
#' @return A dataframe.
#' @import vcfR
#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors unstrsplit
#' @export
bacsnp.transformation <- function(vcf_data){

    fixed_info <- as.data.frame(vcf_data@fix)
    genotype_info <- as.data.frame(vcf_data@gt)

    noISO <- dim(genotype_info)[2]-1
    fixed_info <- do.call(rbind, replicate(noISO, fixed_info, simplify = FALSE))

    split_genotype_info <- function(df) {
      # Extrahieren der Variablennamen aus der FORMAT Spalte
      variable_names <- unlist(str_split(as.character(df$FORMAT[1]), ":"))
      # Erstellen eines leeren Dataframes zum Speichern der Ergebnisse
      result_df <- data.frame(matrix(ncol = length(variable_names) + 1, nrow = 0))
      colnames(result_df) <- c(variable_names, "ISO")

      # Durchlaufen aller Spalten außer 'FORMAT'
      for (col in names(df)[-1]) {
        # Aufteilen der Daten in der Spalte
        split_data <- strsplit(as.character(df[[col]]), ":")
        # Umwandlung der Liste in einen Dataframe
        temp_df <- do.call(rbind, lapply(split_data, function(x) as.data.frame(matrix(x, ncol = length(variable_names), byrow = TRUE))))
        names(temp_df) <- variable_names
        # Hinzufügen der ISO-Spalte
        temp_df$ISO <- col
        # Zusammenführen mit dem Ergebnis-Dataframe
        result_df <- rbind(result_df, temp_df)
      }
      return(result_df)
    }

    genotype_info <- split_genotype_info(genotype_info)

    vcf_info <- cbind(fixed_info, genotype_info)

    vcf_info$POS <- as.numeric(vcf_info$POS)
    vcf_info$COV <- as.numeric(vcf_info$DP)

    if(any(grepl("DPR", colnames(vcf_info))) == TRUE){
      #if DPR is present then do...
      vcf_info$COV.REF <- as.numeric(sapply(strsplit(vcf_info$DPR, ","), function(x) x[1]))
      vcf_info$COV.ALT1 <- vcf_info$COV - vcf_info$COV.REF
      vcf_info$REL.REF <- vcf_info$COV.REF/vcf_info$COV
      vcf_info$REL.ALT1 <- 1-vcf_info$REL.REF

    }else{
      #if DPR is NOT present do...
      vcf_info$COV.REF <- as.numeric(vcf_info$RD)
      vcf_info$COV.ALT1 <- as.numeric(vcf_info$AD)
      vcf_info$REL.REF <- as.numeric(vcf_info$RD)/vcf_info$COV
      vcf_info$REL.ALT1 <- as.numeric(vcf_info$AD)/vcf_info$COV

    }

    return(vcf_info)
  }
