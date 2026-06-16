#' Coupling Coordination Degree (CCD)
#'
#' @param data A numeric matrix or data.frame. Rows are observations,
#' columns are indicators.
#' @param weight Numeric vector of indicator weights. Must have length equal
#' to `ncol(data)`. If `NULL`, equal weights are used.
#' @param method Coupling model. One of `"standard"`, `"wang"`, or `"fan"`.
#' @param threads Number of threads used in computation.
#'
#' @return A data.frame with:
#' \itemize{
#'   \item \code{C}: coupling degree
#'   \item \code{D}: coordination degree
#' }
#' 
#' @export
#'
#' @details
#' The coordination degree is defined as:
#' \deqn{D = \sqrt{C \times T}}
#' where \eqn{T} is the weighted sum of indicators.
#'
#' Full model definitions and formulas are available at:
#' \url{https://github.com/stscl/coupling/discussions/3}
#' 
#' @note
#' Input values should be normalized to `[0, 1]`.
#'
#' @examples
#' set.seed(42)
#' mat = matrix(runif(20), nrow = 5)
#' coupling::ccd(mat)
#'
ccd = \(data, weight = NULL, method = c("standard", "wang", "fan"), threads = 1){
  mat = as.matrix(data)
  method = match.arg(method)
  if (is.null(weight)) weight = rep(1, times = ncol(mat)) / ncol(mat)
  return(RcppCCD(mat, weight, method, threads))
}
