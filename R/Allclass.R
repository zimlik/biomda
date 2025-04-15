#' @title BioMDA class
#' @docType class
#' @slot metabolite The synonyms of compound.
#' @slot metaboliteInfo A list containing biological prior knowledge of query
#' metabolite.
#' @slot disease The synonyms of disease.
#' @slot diseaseInfo A list containing biological prior knowledge of query
#' disease.
#' @slot branch The branch of BioMDA pipeline, one of "tgtp", "ssimm" and
#' "coabm".
#' @slot score Threshold of significance to include a interaction, a number
#' between 0 and 1000 (Highest >= 900, High >= 700, medium >= 400 and low >=
#' 150. Default: 700).
#' @slot node Network nodes represent proteins or metabolites.
#' @slot edge Edges represent protein-protein, compound-protein, or
#' compound-compoud associations.
#' @slot diseaseGeneORA A list containing the result of disease-related genes
#' ORA.
#' @slot enrichmentResult A list of \code{enrichResult} object.
#' imported from DOSE package.
#' @importFrom methods setClass
#' @exportClass BioMDA
setClass("BioMDA",
  slots = c(
    metabolite = "character",
    metaboliteInfo = "list",
    disease = "character",
    diseaseInfo = "list",
    branch = "character",
    score = "numeric",
    node = "data.frame",
    edge = "data.frame",
    diseaseGeneORA = "list",
    enrichmentResult = "list"
  )
)

#' @title BioMDA object
#' @description Construct a BioMDA object.
#' @param metabolite The synonyms of metabolite/compound (ignored case), e.g.
#' "aspirin" or some external id.
#' @param tgtp_score Threshold of significance to search known metabilite-
#' protein interactions (i.e., target proteins of query metabilite), a number
#' between 0 and 1000 (Highest >= 900, High >= 700, medium >= 400 and low >=
#' 150. Default: 700).
#' @param coab_metabolites A list of co-abundant metabolites of query metabolite
#' (i.e., consistency abundance module), identified by WGCNA.
#' @param disease The synonyms of disease, e.g. "Crohn' disease" or "DOID:8778"
#' (Disease Ontology ID). When the exact disease name is uncertain, it may be
#' identified using \code{disease2doid("crohn", fixed = FASLE)} for fuzzy
#' matching.
#' @importFrom methods new
#' @details
#' Some external id of "aspirin":
#' * \code{KEGG compound} D00109
#' * \code{ChEBI} CHEBI:15365
#' * \code{HMDB} HMDB0001879
#' * \code{CAS-RN} 50-78-2
#' * etc.
#'
#' The descriptions of metaboliteInfo:
#' * \code{$cid} the PubChem Cid of metabolite.
#' * \code{$target_proteins} the target proteins of metabolite, accessed by
#' \code{tgtp_search}.
#' * \code{$ssim_metabolites} a set of metabolites shared similar structure with
#' query metabolite, accessed by \code{ssimcid_search}.
#' * \code{$coab_metabolites} a set of co-abundant metabolites of query
#' metabolite (i.e., consistency abundance module), identified by WGCNA.
#'
#' The descriptions of diseaseInfo:
#' * \code{$doid} Disease Ontology ID of disease.
#' * \code{superclass} superclass of disease.
#' * \code{disease_genes} a data.frame of disease-related genes, accessed by
#' \code{disgene_search}.
#'
#' The descriptions of diseaseGeneORA:
#' * \code{$disease} Synonyms of disease.
#' * \code{$geneRatio} (Number of interaction genes in the disease-related
#' genes) / (Total number of interaction genes)
#' * \code{$bgRatio} (Number of disease-related genes)
#' / (Total number of background genes, i.e., universe)
#' * \code{$pvalue} P value of hypergeometric test.
#' * \code{$geneID} The Overlap genes between genes and disease-related genes.
#' * \code{$count} Number of interaction genes in the disease-related genes.
#'
#' The descriptions of enrichmentResult, a list containing:
#' * \code{$`Annotated Keywords (UniProt)`}
#' * \code{$`Biological Process (Gene Ontology)`}
#' * \code{$`Cellular Component (Gene Ontology)`}
#' * \code{$`Molecular Function (Gene Ontology)`}
#' * \code{$`Disease-gene associations (DISEASES)`}
#' * \code{$`Human Phenotype (Monarch)`}
#' * \code{$`KEGG Pathways`}
#' * \code{$`Local Network Cluster (STRING)`}
#' * \code{$`Protein Domains (Pfam)`}
#' * \code{$`Protein Domains (SMART)`}
#' * \code{$`Protein Domains and Features (InterPro)`}
#' * \code{$`Reactome Pathways`}
#' * \code{$`Subcellular localization (COMPARTMENTS)`}
#' * \code{$`Tissue expression (TISSUES)`}
#' * \code{$`WikiPathways`}
#' * \code{$`Disease Gene Network (DisGeNET)`}
#' * \code{$`Network of Cancer Genes (NCG)`}
#' * \code{$`Disease Ontology (DO)`}
#' @return BioMDA object.
#' @seealso \code{BioMDA-class}
#' @examples
#' library(biomda)
#' # Query for diseases containing "crohn".
#' disease2doid("crohn", fixed = FALSE)
#' ##   doid      disease         superclass                 gene_count
#' ##   <chr>     <chr>           <chr>                           <int>
#' ## 1 DOID:8778 Crohn's disease inflammatory bowel disease        110
#' biomda <- BioMDA(metabolite = "Urobilin", disease = "Crohn's disease")
#' @export
BioMDA <- function(metabolite,
                   tgtp_score = 700,
                   coab_metabolites = NULL,
                   disease = NULL) {
  cid <- compound2cid(metabolite)$cid
  tgtp <- tgtp_search(cid, score = tgtp_score)
  ssimcid <- ssimcid_search(cid)
  metaboliteInfo <- list(cid = cid,
                         target_proteins = tgtp,
                         ssim_metabolites = ssimcid)
  if (is.null(coab_metabolites)) {
    char0 <- character(0)
    metaboliteInfo$coab_metabolites <- tibble(coabcid = char0, module = char0,
                                              name = char0, description = char0)
  } else {
    metaboliteInfo$coab_metabolite <- coab_metabolites
  }
  biomda <- new("BioMDA",
                metabolite = metabolite,
                metaboliteInfo = metaboliteInfo)
  if (!is.null(coab_metabolites)) {
    biomda@coab_metabolites <- coab_metabolites
  }
  if (!is.null(disease)) {
    d2doid <- disease2doid(disease)
    diseaseInfo <- list(doid = d2doid$doid,
                        superclass = d2doid$superclass,
                        disease_genes = disgene_search(d2doid$doid))
    biomda@disease <- disease
    biomda@diseaseInfo <- diseaseInfo
  } else {
    biomda@disease <- NA_character_
    biomda@diseaseInfo <- init_diseaseInfo()
  }
  return(biomda)
}

#' @keywords internal
init_metaboliteInfo <- function() {
  char0 <- character(0)
  int0 <- integer(0)
  target_proteins <- tibble(cid = char0, stringid = char0, score = int0,
                            name = char0, entrezid = char0, description = char0)
  ssim_metabolites <- tibble(cid = char0, ssimcid = char0, evidence = char0,
                             name = char0, description = char0)
  coab_metabolites <- tibble(coabcid = char0, module = char0, name = char0,
                             description = char0)
  res <- list(
    target_proteins = target_proteins,
    ssim_metabolites = ssim_metabolites,
    coab_metabolites = coab_metabolites
  )
  return(res)
}

#' @keywords internal
init_diseaseInfo <- function() {
  res <- list(
    doid = NA_character_,
    superclass = NA_character_,
    disease_genes = tibble(ENTREZID = character(0), SYMBOL = character(0))
  )
  return(res)
}
