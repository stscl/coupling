#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include "coupling.h"

// Wrapper function to preform coupling coordination degree analysis
// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame RcppCCD(const Rcpp::NumericMatrix& mat,
                        const Rcpp::NumericVector& weight,
                        const std::string& method = "standard")
{   
    std::vector<std::vector<double>> mat_std = coupling::convert::mat_r2std(mat, true);
    std::vector<double> wt_std = Rcpp::as<std::vector<double>>(weight);

    // Call ccd function
    std::vector<std::vector<double>> res = coupling::ccd::ccd(
        mat_std, wt_std, method
    );
    
    return Rcpp::DataFrame::create( 
        Rcpp::Named("C") = res[0], 
        Rcpp::Named("D") = res[1]
    );
}
