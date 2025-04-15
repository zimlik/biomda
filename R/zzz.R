#' @importFrom utils packageDescription
#' @importFrom pillar style_subtle
.onAttach <- function(libname, pkgname) {
  pkg_version <- packageDescription(pkgname, fields = "Version")
  msg <- paste0(
    pkgname, " v", pkg_version, "  ",
    "For help: https://github.com/zimlik/biomda/issues", "\n\n"
  )
  citation <- paste0(
    "If you use ", pkgname,
    " in published research, please cite the paper:\n\n",
    biomda_citations()
  )
  res <- paste0(msg, citation, suppressmsg(pkgname))
  packageStartupMessage(
    paste0(strwrap(pillar::style_subtle(res)), collapse = "\n")
  )
  return(NULL)
}

biomda_citations <- function() {
  paste(
    "Zheng, Huimin et al. \"In silico method to maximise the biological",
    "potential of understudied metabolomic biomarkers: a study in",
    "pre-eclampsia.\"Gut vol. 73,2 383-385. 5 Jan. 2024,",
    "doi:10.1136/gutjnl-2022-329312.\n\n"
  )
}

suppressmsg <- function(pkgname) {
  paste0(
    "This message can be suppressed by:\n",
    "suppressPackageStartupMessages(library(", pkgname, "))"
  )
}