
#' @title APIs of BioMDA Database
#' @description Submit an HTML form to APIs of BioMDA Database using the POST
#' method.
#' @param resource The resource locator, e.g. "/c2cid/".
#' @param .params A list includes the name-value parameters will process the
#' form.
#' @importFrom jsonlite fromJSON
#' @importFrom RCurl postForm
#' @examples
#' library(biomda)
#' protein <- c("PTCH1", "TP53", "BRCA1", "BRCA2")
#' resp <- biomda_db_api(resource = "/protein2stringid/",
#'                       .params = list(protein = protein))
#' resp
#' ## $p2string
#' ##   protein       stringv12
#' ## 1    TP53 ENSP00000269305
#' ## 2   PTCH1 ENSP00000332353
#' ## 3   BRCA2 ENSP00000369497
#' ## 4   BRCA1 ENSP00000418960
#' @return The text from the HTTP response or NULL.
#' @export
biomda_db_api <- function(resource,
                          .params = list()) {
  host <- getOption("biomda_db_host", "https://www.medam.top")
  resp <- tryCatch(
    postForm(uri = paste0(host, resource), .params = .params),
    error = function(e) NULL
  )
  if (!is.null(resp)) {
    resp <- fromJSON(resp)
  }
  return(resp)
}
