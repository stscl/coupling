## usethis namespace: start
#' @useDynLib coupling, .registration = TRUE
## usethis namespace: end
NULL

.onLoad = \(...) {
  requireNamespace("Rcpp", quietly = TRUE)
}
