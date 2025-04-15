#' @title Show method extensions to show for BioMDA object.
#' @rdname show-methods
#' @param object BioMDA object
#' @importFrom methods show
#' @return Print info of BioMDA object.
#' @export
setMethod("show", signature("object" = "BioMDA"),
  function(object) {
    if (!is.na(object@disease)) {
      disease <- object@disease
    } else {
      disease <- "[Not defined]"
    }
    cat("## The Biomedical Metabolite-Disease Association:", object@metabolite,
        " ~ ", disease, "\n\n")
    cat("Metabolite: ", object@metabolite, "\n")
    cat("|-PubChem CID: ", object@metaboliteInfo$cid, "\n")
    cat("|-Target proteins: ", nrow(object@metaboliteInfo$target_proteins),
        "\n")
    cat("|-Structurally similar metabolites: ",
        nrow(object@metaboliteInfo$ssim_metabolites), "\n")
    cat("|-Co-abundant metabolites: ", length(get_coab_metabolites(object)),
        "\n\n")
    if (!is.na(object@disease)) {
      cat("Disease: ", object@disease, "\n")
      cat("|-Disease Ontology ID: ", object@diseaseInfo$doid, "\n")
      cat("|-Superclass of disease: ", object@diseaseInfo$superclass, "\n")
      cat("|-Counts of disease related genes: ",
          nrow(object@diseaseInfo$disease_genes), "\n\n")
    }
    if (length(object@branch) > 0) {
      cat("Branch: ", object@branch, "\n")
      cat("|-Score: ", object@score, "\n")
      cat("|-Counts of potential interaction proteins: ",
          nrow(filter(object@node, type == "protein")), "\n\n")
    }
  }
)
