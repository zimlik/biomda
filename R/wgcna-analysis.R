#' @title WGCNA modules
#' @description Identify abundance-correlated metabolites (i.e., consistency
#' abundance module) using weighted correlation network analysis (WGCNA).
#' @param abundance Abundance table, a matrix or data frame in which columns are
#' metabolites and rows are samples.
#' @param networkType Network type, one of "signed", "unsigned" and
#' "signed hybrid".
#' @param corMethod Correlation algorithm for network construction, e.g.,
#' "spearman" or "pearson"
#' @param minModuleSize Minimum size of module or cluster, default is 5.
#' @param cutHeight Maximum dissimilarity (i.e., 1-correlation) that qualifies
#' modules for merging.
#' @param nthreads Default 10.
#' @importFrom WGCNA pickSoftThreshold adjacency mergeCloseModules
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom flashClust flashClust
#' @importFrom stats as.dist
#' @return A data.frame with fields:
#' * \code{coabcid} PubChem CID of metabolites in abundance table
#' * \code{module} Consistency abundance module
#' * \code{name} Name of metabolites in abundance table
#' * \code{description} Description of metabolites in abundance table
#' @examples
#' library(biomda)
#' data(IBD.metabolomics)
#' ## Discovery cohort (PRISM)
#' ibd.discovery <- IBD.metabolomics$discovery
#' abunda <- ibd.discovery$abundance
#' wgcna_mods <- wgcna_analysis(abunda)
#' @export
wgcna_analysis <- function(abundance,
                           networkType = "signed",
                           corMethod = "spearman",
                           minModuleSize = 5,
                           cutHeight = 0.25,
                           nthreads = 10) {
  powers <- c(1:10, seq(12, 30, 2))
  sft <- pickSoftThreshold(abundance, powerVector = powers,
                           networkType = networkType) |>
    suppressPrintAndCat()
  softPower <- sft$powerEstimate
  adjMat <- adjacency(datExpr = abundance,
                      power = softPower,
                      type = networkType,
                      corFnc = "cor",
                      corOptions = list(method = corMethod))
  TOM <- TOMsimilarity_parallel(adjMat, nthreads)
  dissTOM <- 1 - TOM
  metaboTree <- flashClust(as.dist(dissTOM), method = "complete")
  moduleLabels <- cutreeDynamic(dendro = metaboTree,
                                distM = dissTOM,
                                method = "hybrid",
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
  merge <- mergeCloseModules(
    exprData = abundance,
    colors =  moduleLabels,
    corFnc = "cor",
    corOptions = list(method = corMethod),
    cutHeight = cutHeight
  )
  metabolite <- colnames(abundance)
  c2cid <- compound2cid(metabolite) |>
    rename(name = compound, coabcid = cid)
  coabcid <- na.omit(c2cid$coabcid) |>
    sub("^CID", "", x = _)
  mod <- as.character(merge$colors)
  res <- tibble(module = mod, name = metabolite) |>
    left_join(c2cid, by = "name") |>
    relocate(coabcid, .before = module)
  resp <- biomda_db_api(resource = "/cid2info/", .params = list(cid = coabcid))
  if (!is.null(resp)) {
    cnodes <- select(resp$cnodes, id, description)
    res <- res |>
      left_join(cnodes, by = c("coabcid" = "id"))
  } else {
    res <- mutate(res, description = NA_character_)
  }
  return(res)
}