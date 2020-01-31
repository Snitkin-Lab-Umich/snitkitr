% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse.R
\name{get_anc_alleles}
\alias{get_anc_alleles}
\title{Get ancestral state of alleles}
\usage{
get_anc_alleles(tree, mat)
}
\arguments{
\item{tree}{rooted tree}

\item{mat}{allele matrix (rows are variants, columns are samples)}
}
\value{
matrix of most likely ancestral allele for each row in allele matrix and probability that that is the ancestral state
}
\description{
Rereference alleles based on rooted tree
}