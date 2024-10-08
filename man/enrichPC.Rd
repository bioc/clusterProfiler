% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathwayCommons.R
\name{enrichPC}
\alias{enrichPC}
\title{enrichPC}
\usage{
enrichPC(gene, ...)
}
\arguments{
\item{gene}{a vector of genes (either hgnc symbols or uniprot IDs)}

\item{...}{additional parameters, see also the parameters supported by the enricher() function}
}
\value{
A \code{enrichResult} instance
}
\description{
ORA analysis for Pathway Commons
}
\details{
This function performs over-representation analysis using  Pathway Commons
}
