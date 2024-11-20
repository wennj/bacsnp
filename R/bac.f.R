#' Example dataset for bacsnp
#'
#' A dataset containing SNP frequency data for demonstration purposes. This dataset
#' represents a small subset of SNPs with coverage and alternative allele frequency
#' information.
#'
#' @format A data frame with 6 rows and 19 variables:
#' \describe{
#'   \item{CHROM}{Chromosome or reference identifier (character).}
#'   \item{POS}{Position of the SNP on the chromosome (integer).}
#'   \item{ID}{Unused identifier; typically NA (character).}
#'   \item{REF}{Reference nucleotide (character).}
#'   \item{ALT}{Alternative nucleotide(s), separated by commas (character).}
#'   \item{QUAL}{Quality score of the SNP (character).}
#'   \item{FILTER}{Filter status; typically NA (character).}
#'   \item{INFO}{Additional information about the SNP, such as coverage and allele counts (character).}
#'   \item{GT}{Genotype of the SNP (character).}
#'   \item{PL}{Phred-scaled likelihoods for genotypes (character).}
#'   \item{DP}{Total coverage at the SNP position (character).}
#'   \item{DPR}{Per-sample depth for each allele, separated by commas (character).}
#'   \item{GQ}{Genotype quality (character).}
#'   \item{ISO}{Isolate identifier (character).}
#'   \item{COV}{Total coverage at the SNP position (numeric).}
#'   \item{COV.REF}{Coverage of the reference nucleotide (numeric).}
#'   \item{COV.ALT1}{Coverage of the first alternative nucleotide (numeric).}
#'   \item{REL.REF}{Relative frequency of the reference nucleotide (numeric).}
#'   \item{REL.ALT1}{Relative frequency of the first alternative nucleotide (numeric).}
#' }
#'
#' @examples
#' data(bac.f)
#' head(bac.f)
"bac.f"
