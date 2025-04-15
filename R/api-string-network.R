#' @title STRING network
#' @description STRING proterin-protein interactions network.
#' @param stringid STRING id, e.g. ENSP00000269305.
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (Highest >= 900, High >= 700, medium >= 400 and low >=
#' 150. Default: 700).
#' @importFrom tibble add_column
#' @details The synonyms of protein can be converted to STRING id using
#' \code{\link{protein2stringid}}.
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
#' # Use the target proteins branch
#' # The ononetin and its catabolites inhibit TRPM3 and activate TRPA1.
#' stringid <- protein2stringid(c("TRPM3", "TRPA1"))$stringid
#' stringnw <- string_network(stringid, score = 700)
#' stringnw$edges
#' stringnw$nodes
#' @export
string_network <- function(stringid, score = 700) {
  resp <- biomda_db_api(resource = "/stringppi/v12/",
                        .params = list(stringid = stringid, score = score))
  if (!is.null(resp)) {
    edges <- resp$edges |>
      mutate(type = "ppi", structure = 0) |>
      select(node1, node2, type, score, experiment, database, neighborhood,
             fusion, cooccurence, structure, coexpression, textmining)
    nodes <- resp$nodes
    network <- list(edges = edges, nodes = nodes)
  } else {
    network <- list(edges = NULL, nodes = NULL)
  }
  return(network)
}
