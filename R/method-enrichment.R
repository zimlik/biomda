#' @title Enrichment analysis
#' @description Enrichment analysis for the BioMDA object using protein
#' associated terms obtained from STRING database (version 12).
#' @param obj BioMDA object containing result of branch, see
#' \code{\link{biomda_branch}}.
#' @param padj P-values adjusted methods, one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param pcutoff 0.05.
#' @importFrom clusterProfiler enricher
#' @references Szklarczyk, Damian et al. “The STRING database in 2023:
#' protein-protein association networks and functional enrichment analyses for
#' any sequenced genome of interest.” Nucleic acids research vol. 51,D1 (2023):
#' D638-D646. doi:10.1093/nar/gkac1000.
#' @examples
#' library(biomda)
#' biomda <- BioMDA("Urobilin", "Crohn's disease")
#' ssimm_biomda <- biomda_branch(biomda, branch = "ssimm", score = 700)
#' ssimm_biomda <- biomda_enrich(ssimm_biomda)
#' @return BioMDA object containing result of enrichment analysis.
#' @details A list of \code{enrichResult} object.
#' @export
setGeneric("biomda_enrich",
  function(obj, padj = "BH", pcutoff = 0.05) standardGeneric("biomda_enrich"),
  signature = "obj"
)

#' @aliases biomda_enrich,BioMDA
#' @rdname biomda_enrich
#' @export
setMethod("biomda_enrich", signature("obj" = "BioMDA"),
  function(obj, padj = "BH", pcutoff = 0.05) {
    if (!inherits(obj@node, "data.frame")) {
      error_msg <- paste0(
        "The slot '@node' of 'BioMDA' object is NULL!\n\n",
        "Use 'biomda_branch(obj, branch = c(\"tgtp\", \"ssimm\", \"coabm\"), ",
        "score = 700)' to get potential proteins of metabolite."
      )
      stop(error_msg)
    }
    biomda_terms <- get_biomda_terms()
    pnode <- obj@node |>
      filter(type == "protein")
    stringid <- pull(pnode, id)
    entrezid <- pull(pnode, entrezid)
    eg_terms_category <- c("KEGG Pathways",
                           "Disease Ontology",
                           "Disease Ontology (DO)",
                           "Network of Cancer Genes (NCG)",
                           "Disease Gene Network (DisGeNET)")
    enres_list <- lapply(names(biomda_terms), function(category) {
      if (category %in% eg_terms_category) {
        genelist <- entrezid
        genetype <- "ENTREZID"
      } else {
        genelist <- stringid
        genetype <- "ENSP"
      }
      gsid2gene <- biomda_terms[[category]]$gsid2gene
      gsid2name <- biomda_terms[[category]]$gsid2name
      enres <- enricher(genelist, pAdjustMethod = padj, pvalueCutoff = pcutoff,
                        TERM2GENE = gsid2gene, TERM2NAME = gsid2name)
      enres <- geneid2name(enres, pnode, genetype)
    })
    names(enres_list) <- names(biomda_terms)
    obj@enrichmentResult <- enres_list
    return(obj)
  }
)