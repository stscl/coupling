#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include "coupling.h"

// Wrapper function to preform coupling coordination degree analysis
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppCCD(const Rcpp::NumericMatrix& mat,
                   const Rcpp::NumericVector& weight,
                   const std::string& method = "standard")
{   
    std::vector<std::vector<double>> mat_std = coupling::convert::mat_r2std(mat, true);
    std::vector<double> wt_std = Rcpp::as<std::vector<double>>(weight);
    
}