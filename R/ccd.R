ccd = \(data, weight = NULL, method = c("standard", "wang", "fan"), threads = 1){
  mat = as.matrix(data)
  method = match.arg(method)
  if (is.null(weight)) weight = rep(1, times = ncol(mat)) / ncol(mat)
  return(RcppCCD(mat, weight, method, threads))
}
