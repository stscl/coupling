metacoupling = \(data, swm_peri = NULL, swm_tele = NULL, weight = NULL, 
                 method = c("standard", "wang", "fan"), threads = 1){
  mat = as.matrix(data)
  method = match.arg(method)
  if (is.null(weight)) weight = rep(1, times = ncol(mat)) / ncol(mat)
  if (is.null(swm_peri)) swm_peri = matrix(0, nrow(data), nrow(data))
  if (is.null(swm_tele)) swm_tele = matrix(0, nrow(data), nrow(data))
  return(RcppMetaCoupling(mat, swm_peri, swm_tele, weight, method, threads))
}
