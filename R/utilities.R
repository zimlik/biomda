## Convert string id or entrez gene id to gene symbol for \code{enrichResult}.
#' @keywords internal
geneid2name <- function(enrichres, pnode, id = c("ENSP", "ENTREZID")) {
  if (inherits(enrichres, "enrichResult")) {
    res <- enrichres@result
    geneid <- strsplit(res$geneID, split = "/")
    gename <- lapply(geneid, function(x) {
      if (id == "ENSP") {
        idx <- match(x, pnode$id)
      } else {
        idx <- match(x, pnode$entrezid)
      }
      paste0(pnode$name[idx], collapse = "/")
    })
    res$geneID <- unlist(gename)
    enrichres@result <- res
  }
  return(enrichres)
}

## Convert synonyms of compound to md5sum string.
#' @importFrom digest digest
#' @importFrom tibble tibble
#' @keywords internal
compound2md5sum <- function(compound) {
  md5sum <- compound |>
    tolower() |>
    lapply(digest, algo = "md5", serialize = FALSE) |>
    unlist()
  compound2md5sum <- tibble(compound = compound, md5sum = md5sum)
  return(compound2md5sum)
}

## Get a set of co-abundant metabolites of query metabolite.
#' @keywords internal
get_coab_metabolites <- function(obj) {
  metabolite <- obj@metabolite
  coabres <- obj@metaboliteInfo$coab_metabolites
  if (nrow(coabres) == 0) {
    coabcid <- character(0)
  } else {
    mod <- coabres |>
      filter(name == metabolite) |>
      pull(module)
    coabcid <- coabres |>
      filter(module == mod) |>
      pull(coabcid)
  }
  return(coabcid)
}