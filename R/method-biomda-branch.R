#' @title BioMDA branch
#' @description The pipeline of BioMDA to predict potential metabolite-disease
#' associations in three different branches. Ref: fig 1.C in Huimin Zheng, et
#' al. (2022).
#' @importFrom methods setGeneric
#' @param obj BioMDA object
#' @param branch The branch name, one of "tgtp", "ssimm" and "coabm".
#' @param score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (Highest >= 900, High >= 700, medium >= 400 and low >=
#' 150. Default: 700).
#' @details The branch of BioMDA pipeline:
#' * \code{tgtp} the target proteins branch.
#' * \code{ssimm} the structurally similar metabolites branch.
#' * \code{coabm} the co-abundant metabolites branch.
#' @references Huimin Zheng, et al. (2022). In silico method to maximise the
#' biological potential of understudied metabolomic biomarkers: a study in
#' pre-eclampsia.
#' @examples
#' library(biomda)
#' biomda <- BioMDA("Urobilin", disease == "Crohn's disease")
#' ssimm_biomda <- biomda_branch(biomda, branch = "ssimm", score = 700)
#' @return BioMDA object containing result of branch.
#' @export
setGeneric("biomda_branch",
  function(obj, branch = c("tgtp", "ssimm", "coabm"), score = 700) {
    standardGeneric("biomda_branch")
  },
  signature = "obj"
)

#' @aliases biomda_branch,BioMDA
#' @rdname biomda_branch
#' @importFrom methods setMethod
#' @export
setMethod("biomda_branch", signature("obj" = "BioMDA"),
  function(obj, branch = c("tgtp", "ssimm", "coabm"), score = 700L) {
    branch <- match.arg(branch, c("tgtp", "ssimm", "coabm"))
    if (branch == "tgtp") {
      stringid <- obj@metaboliteInfo$target_proteins |>
        pull(stringid)
      branch <- "the target proteins branch"
      nw <- string_network(stringid = stringid, score = score)
    } else {
      if (branch == "ssimm") {
        cid <- obj@metaboliteInfo$ssim_metabolites |>
          pull(ssimcid)
        branch <- "the structurally similar metabolites branch"
      } else {
        cid <- get_coab_metabolites(obj)
        branch <- "the co-abundant metabolites branch"
      }
      nw <- stitch_network(cid = cid, score = score)
    }
    obj@branch <- branch
    obj@score <- score
    obj@node <- nw$node
    obj@edge <- nw$edge
    return(obj)
  }
)