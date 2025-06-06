% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wgcna-analysis.R
\name{wgcna_analysis}
\alias{wgcna_analysis}
\title{WGCNA modules}
\usage{
wgcna_analysis(
  abundance,
  networkType = "signed",
  corMethod = "spearman",
  minModuleSize = 5,
  cutHeight = 0.25,
  nthreads = 10
)
}
\arguments{
\item{abundance}{Abundance table, a matrix or data frame in which columns are
metabolites and rows are samples.}

\item{networkType}{Network type, one of "signed", "unsigned" and
"signed hybrid".}

\item{corMethod}{Correlation algorithm for network construction, e.g.,
"spearman" or "pearson"}

\item{minModuleSize}{Minimum size of module or cluster, default is 5.}

\item{cutHeight}{Maximum dissimilarity (i.e., 1-correlation) that qualifies
modules for merging.}

\item{nthreads}{Default 10.}
}
\value{
A data.frame with fields:
\itemize{
\item \code{coabcid} PubChem CID of metabolites in abundance table
\item \code{module} Consistency abundance module
\item \code{name} Name of metabolites in abundance table
\item \code{description} Description of metabolites in abundance table
}
}
\description{
Identify abundance-correlated metabolites (i.e., consistency
abundance module) using weighted correlation network analysis (WGCNA).
}
\examples{
library(biomda)
data(IBD.metabolomics)
## Discovery cohort (PRISM)
ibd.discovery <- IBD.metabolomics$discovery
abunda <- ibd.discovery$abundance
wgcna_mods <- wgcna_analysis(abunda)
}
