#' @title Target proteins of metabolite
#' @description Search known metabolite-protein interactions (experimentally
#' determined or from curated databases) in STITCH.
#' @param cid PubChem cid, e.g. 'CID2244' ('aspirin').
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (Highest >= 900, High >= 700, medium >= 400 and low >=
#' 150. Default: 700).
#' @details The synonyms of compound can be converted to PubChem cid using
#' \code{\link{compound2cid}}.
#' @return A tibble containing the queried Compound (cid) along with
#' comprehensive information on target proteins, including STRING database IDs
#' (stringid), interaction scores (score), protein names (name), Entrez gene IDs
#' (entrezid), and detailed protein descriptions.
#' @examples
#' library(biomda)
#' cid <- compound2cid("N-methylserotonin")$cid # "CID150885"
#' cid2tgtp <- tgtp_search(cid)
#' @export
tgtp_search <- function(cid, score = 700) {
  cid <- sub("^CID", "", cid)
  resp <- biomda_db_api(resource = "/targetcpi/v5/",
                        .params = list(cid = cid, score = score))
  if (!is.null(resp)) {
    cid2tgtp <- resp$tgtp_cpi |>
      mutate(node1 = paste0("CID", node1)) |>
      left_join(resp$tgtp_pnode, by = c("node2" = "id")) |>
      rename(cid = node1, stringid = node2) |>
      as_tibble()
  } else {
    cid2tgtp <- tibble(cid = character(0),
                       stringid = character(0),
                       score = integer(0),
                       name = character(0),
                       entrezid = character(0),
                       description = character(0))
  }
  return(cid2tgtp)
}
