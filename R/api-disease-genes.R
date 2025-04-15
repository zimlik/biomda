#' @title Retrieve disease-related genes.
#' @description Retrieve disease-related genes.
#' @param doid Disease Ontology ID (DOID). Note: the length of "doid" must be 1.
#' @importFrom clusterProfiler bitr
#' @details The aliases of disease can be converted to DOID using
#' \code{\link{disease2doid}}.
#' @return A tibble with ncbi entrez gene id and gene symbol.
#' @examples
#' library(biomda)
#' # Query for pre-eclampsia disease (PE, DOID:10591).
#' pe_disge <- disgene_search("DOID:10591")
#' # Query for morbid obesity disease (DOID:11981).
#' morbid_obesity_disge <- disgene_search("DOID:11981")
#' @export
disgene_search <- function(doid) {
  if (length(doid) > 1) {
    stop("The length of 'doid' must be 1")
  }
  init_disge <- tibble(ENTREZID = NA_character_, SYMBOL = NA_character_)
  resp <- biomda_db_api(resource = "/doid2gene/", .params = list(doid = doid))
  if (!is.null(resp)) {
    eg <- resp$doid2gene |>
      pull(entrezid) |>
      strsplit(split = ",") |>
      unlist()
    disge <- suppressWarnings(
      suppressWarnings(
        bitr(eg, "ENTREZID", "SYMBOL", "org.Hs.eg.db", drop = FALSE) |>
          as_tibble()
      )
    )
  } else {
    disge <- init_disge
  }
  return(disge)
}
