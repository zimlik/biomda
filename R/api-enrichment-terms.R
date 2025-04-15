## Download and cache enrichment terms.
#' @importFrom utils download.file
#' @keywords internal
get_biomda_terms <- function() {
  biomda_env <- get_biomda_env()
  if (!exists("biomda_terms", envir = biomda_env, inherits = FALSE)) {
    url <- "http://39.108.209.221/biomda-enrichment-terms.rds"
    terms_data <- tempfile(fileext = ".rds")
    on.exit(file.remove(terms_data))
    download.file(url, terms_data)
    biomda_terms <- tryCatch(readRDS(terms_data), error = function(e) NULL)
    if (is.null(biomda_terms)) {
      stop("Enrichment terms download failed!")
    } else {
      assign("biomda_terms", biomda_terms, envir = biomda_env)
    }
  }
  biomda_terms <- get("biomda_terms", envir = biomda_env)
  return(biomda_terms)
}

## Create a .biomda_env env to cache enrichment terms.
#' @keywords internal
get_biomda_env <- function() {
  if (!exists(".biomda_env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".biomda_env", new.env(), envir = envir)
  }
  get(".biomda_env", envir = .GlobalEnv)
}
