#' @title Differential abundance analysis
#' @description Differential abundance analysis for metabolomics.
#' @param abundance Abundance table, a matrix or data frame in which columns are
#' metabolites and rows ar samples.
#' @param resp Response to be modelled, a factor (same length as row number of
#' \code{abundance}) with two levels.
#' @param compared Define the condition A and condition B in "response" for
#' calculating fold change. For example, set "compared" to c("case", "control"),
#' then the fold change is "case/control".
#' @param ortho FALSE (for PLS-DA) or TRUE (for OPLS-DA).
#' @param log10L Whether the 'Abudnace Table' table needs to be
#' log10-transformed.
#' @param scaling Scaling methods for "abundance", one of "none", "center",
#' "pareto", and "standard" (default). See details.
#' @param test t.test or wilcox.test.
#' @param padj_method  P-values adjusted methods, one of "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr" and "none".
#' @param vip_cutoff Threshold of VIP score of the (O)PLS-DA model.
#' @param padj_cutoff Threshold of P-value adjusted for t.test (or wilcox.test).
#' @param auc_cutoff Threshold of auc of ROC curve.
#' @param log2fc_cutoff Threshold of absolute value of log2 fold change.
#' @param ... further arguments to be passed to or from methods.
#' @importFrom ropls opls getVipVn
#' @importFrom stats t.test wilcox.test p.adjust na.omit
#' @importFrom pROC roc
#' @importFrom dplyr if_else
#' @details Scaling methods as following:
#' * \code{none} none
#' * \code{center}  mean-centering scaling
#' * \code{pareto} mean-centering and pareto scaling
#' * \code{standard} mean-centering and unit variance scaling, default
#' @return A data.frame with fields:
#' * \code{metabolite} Metabolites
#' * \code{vip} Variable Importance in Projection (VIP) of (O)PLS-DA module
#' * \code{auc} AUC of ROC curve analysis
#' * \code{pval, padj} P value and p value adjusted of t.test (or wilcox.test)
#' * \code{fc, log2fc, fc_direct} Fold change (FC), log2(FC) and direction of
#' FC (e.g., case/control).
#' * \code{significant} 1 for significantly differential metabolites, or 0 for
#' not significant.
#' @note Fold change was calculated using the \code{original abundance} table,
#' while p-values and AUC values were computed based on log10-transformed (when
#' \code{log10L = TRUE}, or FALSE for untransformed) abundance data. VIP scores
#' were derived from scaled data (default: \code{scaling = "standard"}, see
#' details) using either log10-transformed or untransformed abundance values.
#' @examples
#' library(biomda)
#' data(IBD.metabolomics)
#' ## Discovery cohort (PRISM)
#' ibd.discovery <- IBD.metabolomics$discovery
#' abunda <- ibd.discovery$abundance
#' metada <- ibd.discovery$metadata
#' # CD v.s. Control
#' metada <- metada[metada$Diagnosis != "UC", ]
#' table(metada$Diagnosis) # 68 Crohn' disease and 34 control samples
#' abunda <- abunda[match(rownames(metada), rownames(abunda)), ]
#' # OPLS-DA
#' compared <- c("CD", "Control")
#' diffres <- diff_metabolites(abunda, metada$Diagnosis, compared, ortho = TRUE)
#' @export
diff_metabolites <- function(abundance,
                             resp,
                             compared,
                             ortho = FALSE,
                             log10L = FALSE,
                             scaling = "standard",
                             test = "t.test",
                             padj_method = "BH",
                             vip_cutoff = 1,
                             padj_cutoff = 0.05,
                             auc_cutoff = 0.7,
                             log2fc_cutoff = 1.5,
                             ...) {
  lev <- na.omit(unique(resp))
  if (length(lev) != 2) {
    stop("The level of response factor is not equal 2.")
  }
  if (missing(compared) || length(compared) != 2 || !all(compared %in% lev)) {
    compared <- sort(unique(resp))
    patt <- c("The 'compared' is missing or",
              "is not equal to the levels of 'resp',",
              "then set to c(%s) automatically.")
    patt <- paste0(patt, collapse = " ")
    mesg <- paste0("\"", compared, "\"", collapse = ",")
    mesg <- sprintf(patt, mesg)
    warning(mesg)
  }
  if (is.null(list(...)$origin)) {
    origin_abun <- abundance
  } else {
    origin_abun <- list(...)$origin
  }
  if (log10L) {
    abundance <- log10(abundance + 10^-6)
  }
  if (ortho) {
    oplsda <- opls(abundance, resp, predI = 1, orthoI = NA, scaleC = scaling,
                   fig.pdfC = "none") |>
      suppressPrintAndCat()
    vip <- getVipVn(oplsda)
  } else {
    plsda <- opls(abundance, resp, scaleC = scaling, fig.pdfC = "none") |>
      suppressPrintAndCat()
    vip <- getVipVn(plsda)
  }
  statistic <- apply(abundance, 2, function(x) {
    if (test == "t.test") {
      fit <- t.test(x ~ resp)
    } else if (test == "wilcox.test") {
      fit <- wilcox.test(x ~ resp)
    }
    roc_curve <- suppressMessages(roc(resp, x, levels = rev(compared)))
    auc_value <- roc_curve$auc
    list(pval = fit$p.value, auc_value = auc_value)
  })
  pval <- lapply(statistic, `[[`, "pval") |> unlist()
  adj_pval <- p.adjust(pval, method = padj_method)
  auc_value <- lapply(statistic, `[[`, "auc_value") |> unlist()
  fc <- apply(origin_abun, 2, function(x) {
    m <- tapply(x, resp, mean, na.rm = TRUE)
    m[compared[1]] / m[compared[2]]
  })
  log2fc <- log2(fc)
  fc_direct <- paste0(compared, collapse = "/")
  metabolite <- colnames(abundance)
  res <- tibble(metabolite = metabolite, vip = vip, pvalue = pval,
                padj = adj_pval, auc = auc_value, fc = fc, log2fc = log2fc,
                fc_direct = fc_direct)
  auc_cutoff2 <- 1 - auc_cutoff
  res <- res |>
    mutate(significant = as.integer(vip > vip_cutoff &
                                    padj < padj_cutoff &
                                    abs(log2fc) > log2fc_cutoff &
                                    (auc > auc_cutoff | auc < auc_cutoff2))) |>
    mutate(significant = if_else(is.na(significant), 0, significant))

  return(res)
}