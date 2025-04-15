#' @title Disease-related genes ORA
#' @description Determine whether genes from disease-related genes are enriched
#' more than would be expected in potential interaction genes of metabolite.
#' Ref: fig 1.C in Huimin Zheng, et al. (2022).
#' @param obj BioMDA object containing result of branch, see
#' \code{\link{biomda_branch}}.
#' @param universe Count of universal background protein-coding genes.
#' For example, background of homo sapiens is 21306.
#' @importFrom stats phyper
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @return BioMDA object containing result of disease-related ORA.
#' @details A list contained:
#' * \code{disease} Synonyms of disease.
#' * \code{geneRatio} (Number of interaction genes in the disease-related genes)
#' / (Total number of interaction genes)
#' * \code{bgRatio} (Number of disease-related genes)
#' / (Total number of background genes, i.e., universe)
#' * \code{pvalue} P value of hypergeometric test.
#' * \code{geneID} The Overlap genes between genes and disease-related genes.
#' * \code{count} Number of interaction genes in the disease-related genes.
#' @examples
#' library(biomda)
#' biomda <- BioMDA("Urobilin", "Crohn's disease")
#' ssimm_biomda <- biomda_branch(biomda, branch = "ssimm", score = 700)
#' ssimm_biomda <- biomda_disge_ora(ssimm_biomda)
#' @export
setGeneric("biomda_disge_ora",
  function(obj, universe = 21306) standardGeneric("biomda_disge_ora"),
  signature = "obj"
)

#' @aliases biomda_disge_ora,BioMDA
#' @rdname biomda_disge_ora
#' @export
setMethod("biomda_disge_ora", signature("obj" = "BioMDA"),
  function(obj, universe = 21306) {
    if (!inherits(obj@node, "data.frame")) {
      error_msg <- paste0(
        "The slot '@node' of 'BioMDA' object is NULL!\n\n",
        "Use 'biomda_branch(obj, branch = c(\"tgtp\", \"ssimm\", \"coabm\"), ",
        "score = 700)' to get potential proteins of metabolite."
      )
      stop(error_msg)
    }
    node <- obj@node
    disgene <- obj@diseaseInfo$disease_genes |>
      pull(ENTREZID)
    disease <- obj@disease
    bgRatio <- paste0(length(disgene), "/", universe)
    res <- list(disease = disease, geneRatio = "0/0", bgRatio = bgRatio,
                pvalue = NA, geneID = "", count = 0)
    if (!is.null(node)) {
      iprot <- node |> filter(type == "protein")
      entrezid <- pull(iprot, entrezid)
      q <- sum(entrezid %in% disgene) - 1
      m <- length(disgene)
      n <- universe - m
      k <- length(entrezid)
      if (q > 0) {
        res$pvalue <- phyper(q, m, n, k, lower.tail = FALSE)
        geneID <- iprot |>
          filter(entrezid %in% disgene) |>
          pull(name)
        res$count <- length(geneID)
        res$geneID <- paste0(geneID, collapse = "/")
        res$geneRatio <- paste0(q + 1, "/", k)
      }
    }
    obj@diseaseGeneORA <- res
    return(obj)
  }
)
