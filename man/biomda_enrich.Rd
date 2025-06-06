% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method-enrichment.R
\name{biomda_enrich}
\alias{biomda_enrich}
\alias{biomda_enrich,BioMDA-method}
\alias{biomda_enrich,BioMDA}
\title{Enrichment analysis}
\usage{
biomda_enrich(obj, padj = "BH", pcutoff = 0.05)

\S4method{biomda_enrich}{BioMDA}(obj, padj = "BH", pcutoff = 0.05)
}
\arguments{
\item{obj}{BioMDA object containing result of branch, see
\code{\link{biomda_branch}}.}

\item{padj}{P-values adjusted methods, one of "holm", "hochberg", "hommel",
"bonferroni", "BH", "BY", "fdr", "none".}

\item{pcutoff}{0.05.}
}
\value{
BioMDA object containing result of enrichment analysis.
}
\description{
Enrichment analysis for the BioMDA object using protein
associated terms obtained from STRING database (version 12).
}
\details{
A list of \code{enrichResult} object.
}
\examples{
library(biomda)
biomda <- BioMDA("Urobilin", "Crohn's disease")
ssimm_biomda <- biomda_branch(biomda, branch = "ssimm", score = 700)
ssimm_biomda <- biomda_enrich(ssimm_biomda)
}
\references{
Szklarczyk, Damian et al. “The STRING database in 2023:
protein-protein association networks and functional enrichment analyses for
any sequenced genome of interest.” Nucleic acids research vol. 51,D1 (2023):
D638-D646. doi:10.1093/nar/gkac1000.
}
