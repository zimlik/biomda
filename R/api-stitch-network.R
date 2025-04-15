#' @title STITCH network
#' @description STITCH compound-protein, compound-compoud, interactions network.
#' @param cid PubChem cid, e.g. 'CID2244' ('aspirin').
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (Highest >= 900, High >= 700, medium >= 400 and low >=
#' 150. Default: 700).
#' @importFrom dplyr select mutate_if bind_rows
#' @importFrom tidyr replace_na
#' @details The synonyms of compound can be converted to PubChem cid using
#' \code{\link{compound2cid}}.
#' @return A list containing weights of each edge and the annotations of each
#' node as following:
#'
#' the tibble \code{edges} with fields:
#' * \code{node1, node2} protein- (ENSPxxx) or compound-nodes (CIDxxx) in the
#' network.
#' * \code{type} type of interaction edges in network, e.g., cpi
#' (compound-protein interactions), cci (compound-compound interactions) and ppi
#' (protein-protein interactions).
#' * \code{score} combined score.
#' * \code{experiment} experimentally determined (rank: known interactions).
#' * \code{database} from curated databases (rank: known interactions).
#' * \code{neighborhood} gene neighborhood (rank: predicted interactions).
#' * \code{fusion} gene fusions (rank: predicted interactions).
#' * \code{cooccurence} gene co-occurrence (rank: predicted interactions).
#' * \code{structure} extrapolated 'compound-protein associations' based on
#' chemical structure (rank: predicted interactions).
#' * \code{coexpression} co-expression (rank: others).
#' * \code{textmining} textmining (rank: others).
#'
#' the tibble \code{nodes} with fields:
#' * \code{id} PubChem cid (compound) or STRING id (protein).
#' * \code{type} type of nodes in network, protein or metabolite.
#' * \code{name} common synonyms of compound or protein.
#' * \code{external} external id, NA (for compounds) or Entrez Gene ID (for
#' proteins).
#' * \code{description} description of compound or protein.
#' @examples
#' library(biomda)
#' # Use the structurally similar metabolites branch
#' cid <- compound2cid("ononetin") # "CID259632"
#' ssimcid <- ssimcid_search(cid)$ssimcid
#' stitchnw <- stitch_network(ssimcid, score = 700)
#' stitchnw$edges
#' stitchnw$nodes
#' @export
stitch_network <- function(cid, score = 700) {
  cid <- sub("^CID", "", cid)
  resp <- biomda_db_api(resource = "/stitchcpi/v5/",
                        .params = list(cid = cid, score = score))
  if (!is.null(resp)) {
    cpi <- resp$edges$cpi |>
      rename(structure = prediction) |>
      mutate(type = "cpi", node1 = paste0("CID", node1), neighborhood = 0L,
             fusion = 0L, cooccurence = 0L, coexpression = 0L)
    if (inherits(resp$edges$ppi, "data.frame")) {
      ppi <- resp$edges$ppi |>
        mutate(type = "ppi")
    } else {
      ppi <- NULL
    }
    if (inherits(resp$edges$cci, "data.frame")) {
      cci <- resp$edges$cci |>
        mutate(type = "cci", node1 = paste0("CID", node1),
               node2 = paste0("CID", node2))
    } else {
      cci <- NULL
    }
    edges <- bind_rows(cpi, ppi, cci) |>
      mutate_if(is.integer, ~replace_na(., 0)) |>
      select(node1, node2, type, score, experiment, database, neighborhood,
             fusion, cooccurence, structure, coexpression, textmining) |>
      as_tibble()
    pnodes <- resp$nodes$pnodes |>
      mutate(type = "protein")
    cnodes <- resp$nodes$cnodes |>
      mutate(type = "compound", id = paste0("CID", id))
    nodes <- bind_rows(cnodes, pnodes) |>
      select(id, name, type, entrezid, description) |>
      as_tibble()
    network <- list(edges = edges, nodes = nodes)
  } else {
    network <- list(edges = NULL, nodes = NULL)
  }
  return(network)
}
