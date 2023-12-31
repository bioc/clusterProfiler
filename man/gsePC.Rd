% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathwayCommons.R
\name{gsePC}
\alias{gsePC}
\title{gsePC}
\usage{
gsePC(geneList, source, keyType, ...)
}
\arguments{
\item{geneList}{a ranked gene list}

\item{source}{Data source of Pathway Commons, e.g., 'reactome', 'kegg', 'pathbank', 'netpath', 'panther', etc.}

\item{keyType}{specify the type of input 'gene' (one of 'hgnc' or 'uniprot')}

\item{...}{additional parameters, see also the parameters supported by the GSEA() function}
}
\value{
A \code{gseaResult} instance
}
\description{
GSEA analysis for  Pathway Commons
}
\details{
This function performs GSEA using  Pathway Commons
}
