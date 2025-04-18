#' @title Structurally similar compounds
#' @description Search compounds shared similarity structure.
#' @param cid PubChem cid, e.g., CID2244.
#' @importFrom dplyr relocate filter
#' @details \code{ssimcid_search} is based on substructure key-based 2D
#' Tanimoto similarity search (tanimoto >= 90%). The synonyms of compound can be
#' converted to "cid" (PubChem cid) using \code{\link{compound2cid}}.
#' @return A tibble containing the queried cid and ssimcid (compounds shared
#' structural similarity in 2D)
#' @examples
#' library(biomda)
#' cid <- compound2cid("ononetin")$cid # "CID259632"
#' cid2ssimcid <- ssimcid_search(cid)
#' cid2ssimcid$ssimcid
#' @export
ssimcid_search <- function(cid) {
  cid <- sub("^CID", "", cid)
  init_cid2ssimcid <- list(cid = cid,
                           ssimcid = rep(NA_character_, length(cid)))
  resp <- biomda_db_api(resource = "/cid2ssimcid/", .params = list(cid = cid))
  if (!is.null(resp)) {
    cid2ssimcid <- resp$cid2ssimcid |>
      separate_longer_delim(ssimcid, delim = ",") |>
      left_join(resp$cnodes, by = c("ssimcid" = "id")) |>
      filter(!is.na(name), cid != ssimcid) |>
      mutate(cid = paste0("CID", cid),
             ssimcid = paste0("CID", ssimcid),
             evidence = "tanimoto >= 90%") |>
      relocate(evidence, .before = name) |>
      as_tibble()
  } else {
    cid2ssimcid <- as_tibble(init_cid2ssimcid)
  }
  return(cid2ssimcid)
}
