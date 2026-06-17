#' Metacoupling Analysis
#'
#' @inheritParams ccd
#' @param swm_peri A numeric matrix representing the **peri (local) spatial weight matrix**.
#' Must be square with dimension equal to `nrow(data)`. If `NULL`, a zero matrix is used.
#' @param swm_tele A numeric matrix representing the **tele (long-distance) spatial weight matrix**.
#' Must be square with dimension equal to `nrow(data)`. If `NULL`, a zero matrix is used.
#'
#' @return A data.frame with:
#' \itemize{
#'   \item `Intra_C`: intra-system coupling degree
#'   \item `Intra_D`: intra-system coordination degree
#'   \item `Peri_C`: peri-coupling degree
#'   \item `Peri_D`: peri coordination degree
#'   \item `Tele_C`: tele-coupling degree
#'   \item `Tele_D`: tele coordination degree
#' }
#'
#' @export
#'
#' @details
#' Full model definitions and formulas are available at:
#' \url{https://github.com/stscl/coupling/discussions/8}
#'
#' @note
#' Input values should be normalized to `[0, 1]`. Spatial weight matrices are
#' typically symmetric.
#'
#' @examples
#' set.seed(42)
#' mat = matrix(runif(20), nrow = 5)
#' swm1 = matrix(runif(25), 5, 5)
#' swm2 = matrix(runif(25), 5, 5)
#' coupling::metacoupling(mat, swm1, swm2)
#'
metacoupling = \(data, swm_peri = NULL, swm_tele = NULL, weight = NULL, 
                 method = c("standard", "wang", "fan"), threads = 1){
  mat = as.matrix(data)
  method = match.arg(method)
  if (is.null(weight)) weight = rep(1, times = ncol(mat)) / ncol(mat)
  if (is.null(swm_peri)) swm_peri = matrix(0, nrow(data), nrow(data))
  if (is.null(swm_tele)) swm_tele = matrix(0, nrow(data), nrow(data))
  return(RcppMetaCoupling(mat, swm_peri, swm_tele, weight, method, threads))
}
