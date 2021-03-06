% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bobR.R
\name{getYeastUniProtInfo}
\alias{getYeastUniProtInfo}
\title{gets yeast UniProt info from UniProt URL and saves to local file}
\usage{
getYeastUniProtInfo(
  URL = "https://www.uniprot.org/docs/yeast.txt",
  path = "uniprot.txt"
)
}
\arguments{
\item{URL}{The web location of the file.
Defaults to \code{'https://www.uniprot.org/docs/yeast.txt'}}

\item{path}{The directory path to store the file
Defaults to \code{"uniprot.txt"} in working directory.}
}
\description{
This function simply loads a list of the yeast genes in the Swiss UNIPROT database
and saves it as a file
}
