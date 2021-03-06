% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp_filter.R
\name{bacsnp.filter}
\alias{bacsnp.filter}
\title{Filtering SNP frequency dataframe}
\usage{
bacsnp.filter(vcfdf, min.abs.cov = 0, min.abs.alt = 0,
  min.rel.alt = 0)
}
\arguments{
\item{vcfdf}{A dataframe with SNP frequencies. Data was created with vcf_transformation.}

\item{min.abs.cov}{The minimal absolute coverage (COV) that is required to consider SNP positions. All positions with a miniam absolute coverage below this value will be removed.}

\item{min.abs.alt}{An absolute value for the minimal number of occurances of an alternative nucleotide.}

\item{min.rel.alt}{An relative value for the minimal frequency of an alternative nucleotide.}
}
\value{
A vcf dataframe that was filtered by a required minimal total coverage (COV), minimal coverage of alternative nucleotides (COV.ALT) or a minimal relative frequency of the first alterntative (REL.ALT1).
}
\description{
Filtering SNP frequency dataframe
}
\examples{
bacsnp.filter(bac, min.abs.cov = 100)
bacsnp.filter(bac, min.abs.cov = 100, min.abs.alt = 10)
bacsnp.filter(bac, min.abs.cov = 100, min.abs.alt = 10, min.rel.alt = 0.05)
bacsnp.filter(bac, min.rel.alt = 0.05)
}
