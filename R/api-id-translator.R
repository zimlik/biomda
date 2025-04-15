#' @title Convert synonyms of compound to PubChem cid.
#' @description Convert synonyms of compound to PubChem cid.
#' @param compound Synonyms of compound (ignored case), e.g. "aspirin" or
#' some external id (see details).
#' @importFrom tibble as_tibble
#' @importFrom tidyr separate_longer_delim
#' @importFrom dplyr pull mutate group_by arrange slice_head ungroup left_join
#' @details some external id of "aspirin":
#' * \code{KEGG compound} D00109
#' * \code{ChEBI} CHEBI:15365
#' * \code{HMDB} HMDB0001879
#' * \code{CAS-RN} 50-78-2
#' * etc.
#' @examples
#' library(biomda)
#' compound <- c("N-methylserotonin", "Urobilin", "Ononetin")
#' c2cid <- compound2cid(compound = compound)
#' @return A data.frame containing queried compound and PubChem cid.
#' @export
compound2cid <- function(compound) {
  init_c2cid <- list(compound = compound,
                     cid = rep(NA_character_, length(compound)))
  is_cid <- grepl("^CID\\d+$", compound)
  cid <- sub("^CID", "", compound[is_cid])
  if (length(cid) > 0) {
    init_c2cid$cid[is_cid] <- cid
  }
  non_cid <- compound[!is_cid]
  if (length(non_cid) > 0) {
    c2md5sum <- compound2md5sum(non_cid)
    md5sum <- pull(c2md5sum, md5sum)
    resp <- biomda_db_api(resource = "/c2cid/", .params = list(md5sum = md5sum))
    if (!is.null(resp)) {
      resp_c2cid <- resp$md5sum2cid |>
        separate_longer_delim(cid, delim = ",") |>
        group_by(md5sum) |>
        arrange(as.integer(cid), .by_group = TRUE) |>
        slice_head(n = 1) |>
        ungroup() |>
        left_join(c2md5sum, by = "md5sum") |>
        mutate(cid = paste0("CID", cid))
      idx <- match(resp_c2cid$compound, init_c2cid$compound)
      init_c2cid$cid[idx] <- resp_c2cid$cid
    }
  }
  c2cid <- as_tibble(init_c2cid)
  return(c2cid)
}

#' @title Convert synonyms of protein to STRING id.
#' @description Convert synonyms of protein to STRING id.
#' @param protein Synonyms of protein(ignored case), e.g. "TP53".
#' @importFrom dplyr rename
#' @return A data.frame containing queried protein and STRING id.
#' @examples
#' library(biomda)
#' protein <- c("PTCH1", "TP53", "BRCA1", "BRCA2")
#' p2stringid <- protein2stringid(protein = protein)
#' @export
protein2stringid <- function(protein) {
  init_p2string <- list(protein = protein,
                        toupper_protein = toupper(protein),
                        stringid = rep(NA_character_, length(protein)))
  resp <- biomda_db_api(resource = "/protein2stringid/",
                        .params = list(protein = protein))
  if (!is.null(resp)) {
    resp_p2string <- resp$p2string |>
      rename(stringid = id) |>
      mutate(toupper_protein = toupper(protein)) |>
      group_by(toupper_protein) |>
      arrange(stringid, .by_group = TRUE) |>
      slice_head(n = 1) |>
      ungroup()
    idx <- match(resp_p2string$toupper_protein, init_p2string$toupper_protein)
    init_p2string$stringid[idx] <- resp_p2string$stringid
  }
  p2string <- as_tibble(init_p2string[c("protein", "stringid")])
  if (!is.null(resp$pnodes)) {
    p2string <- p2string |>
      left_join(resp$pnodes, by = c("stringid" = "id"))
  }
  return(p2string)
}

#' @title Convert aliases of disease to Diease Ontology id (DOID).
#' @description Convert aliases of disease to Diease Ontology id (DOID).
#' @param disease Aliases of disease(ignored case), e.g. "pre-eclampsia". Note:
#' the length of "disease" must be 1.
#' @param fixed	If TRUE match "disease" exactly, otherwise query "disease" fuzzy
#' faintly.
#' @examples
#' library(biomda)
#' # Query for pre-eclampsia disease (PE).
#' disease2doid("pre-eclampsia")
#' # Query by directly entering the DOID.
#' disease2doid("DOID:10591")
#' # Query for diseases containing "obesity".
#' disease2doid("obesity", fixed = FALSE)
#' @return A data.frame containing DOID, name, superclass and number of related
#' genes of the queried disease.
#' @export
disease2doid <- function(disease, fixed = TRUE) {
  if (length(disease) > 1) {
    stop("The length of 'disease' must be 1")
  }
  init_d2doid <- tibble(doid = NA_character_,
                        disease = disease,
                        superclass = NA_character_,
                        count = 0)
  resp <- biomda_db_api(resource = "/disease2doid/",
                        .params = list(disease = disease, fixed = fixed))
  if (!is.null(resp)) {
    d2doid <- resp$d2doid |>
      rename(count = ngene) |>
      as_tibble()
  } else {
    d2doid <- init_d2doid
  }
  return(d2doid)
}
