#' SNVs after BCF tool filtered MpileUp VCF data.
#'
#' A dataset containing the SNP positions and their nucleotide frequencies
#' of three Illumina sequenced baculoviruses.
#'
#' A dataset containing SNP positions, quality scores, nucleotide frequencies,
#' and relative frequencies for three alternative nucleotides across isolates of baculoviruses.
#' The dataset provides information useful for SNP analysis, including reference and alternative
#' coverage, as well as relative frequencies for each alternative nucleotide.
#'
#' @format A data frame with 6 rows and 14 variables:
#' \describe{
#'   \item{POS}{Integer. Genomic position of the SNP.}
#'   \item{QUAL}{Numeric. Quality score of the SNP.}
#'   \item{ISO}{Factor. Isolate identifier (e.g., \code{iIso1}, \code{iIso2}).}
#'   \item{GQ}{Integer. Genotype quality score.}
#'   \item{REF}{Factor. Reference nucleotide.}
#'   \item{ALT}{Factor. Alternative nucleotides (up to three alternatives).}
#'   \item{COV}{Integer. Total coverage at the SNP position.}
#'   \item{COV.REF}{Integer. Coverage of the reference nucleotide.}
#'   \item{COV.ALT1}{Integer. Coverage of the first alternative nucleotide.}
#'   \item{COV.ALT2}{Integer. Coverage of the second alternative nucleotide.}
#'   \item{COV.ALT3}{Integer. Coverage of the third alternative nucleotide.}
#'   \item{REL.REF}{Numeric. Relative frequency of the reference nucleotide.}
#'   \item{REL.ALT1}{Numeric. Relative frequency of the first alternative nucleotide.}
#'   \item{REL.ALT2}{Numeric. Relative frequency of the second alternative nucleotide.}
#'   \item{REL.ALT3}{Numeric. Relative frequency of the third alternative nucleotide.}
#' }
#'
#' @examples
#' data(bac.df)
#' head(bac.df)
"bac.df"
