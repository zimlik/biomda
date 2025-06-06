% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api-ssim-compound.R
\name{ssimcid_search}
\alias{ssimcid_search}
\title{Structurally similar compounds}
\usage{
ssimcid_search(cid, include_query = TRUE, n)
}
\arguments{
\item{cid}{PubChem cid, e.g., CID2244.}

\item{include_query}{Whether the structurally similar compounds contains the
query cid.}

\item{n}{The number of structurally analogous compounds was analyzed, with
missing values indicating retention of all entries in the dataset. The order
see \code{details}.}
}
\value{
A tibble containing the queried cid and ssimcid (compounds shared
structural similarity in 2D)
}
\description{
Search compounds shared similarity structure.
}
\details{
\code{ssimcid_search} is based on substructure key-based 2D
Tanimoto similarity search (tanimoto >= 90\%). The synonyms of compound can be
converted to "cid" (PubChem cid) using \code{\link{compound2cid}}. The
\code{ssimcid} are sorted by "annotation order", a sort of measure about how
much is known about a compound based on how much.
information we have about it.
}
\examples{
library(biomda)
cid <- compound2cid("ononetin")$cid # "CID259632"
cid2ssimcid <- ssimcid_search(cid)
cid2ssimcid$ssimcid
}
