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
                        const std::string& method = "standard",
                        int threads = 1)
{   
    std::vector<std::vector<double>> mat_std = coupling::convert::mat_r2std(mat, true);
    std::vector<double> wt_std = Rcpp::as<std::vector<double>>(weight);

    // Call ccd function
    std::vector<std::vector<double>> res = coupling::ccd::ccd(
        mat_std, wt_std, method, static_cast<size_t>(std::abs(threads))
    );
    
    return Rcpp::DataFrame::create( 
        Rcpp::Named("C") = res[0], 
        Rcpp::Named("D") = res[1]
    );
}

// Wrapper function to preform meta-coupling analysis
// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame RcppMetaCoupling(
    const Rcpp::NumericMatrix& mat,  
    const Rcpp::NumericMatrix& swm_peri, 
    const Rcpp::NumericMatrix& swm_tele,      
    const Rcpp::NumericVector& weight,                 
    const std::string& method = "standard",                
    int threads = 1)
{   
    std::vector<std::vector<double>> mat_std = coupling::convert::mat_r2std(mat, true);
    std::vector<std::vector<double>> swm_peri_std = coupling::convert::mat_r2std(swm_peri, true);
    std::vector<std::vector<double>> swm_tele_std = coupling::convert::mat_r2std(swm_tele, true);
    std::vector<double> wt_std = Rcpp::as<std::vector<double>>(weight);

    // Call metacoupling function
    std::vector<std::vector<double>> res = coupling::metacoupling::metacoupling(
        mat_std, swm_peri_std, swm_tele_std, wt_std, method, static_cast<size_t>(std::abs(threads))
    );
    
    return Rcpp::DataFrame::create( 
        Rcpp::Named("Intra_C") = res[0], 
        Rcpp::Named("Intra_D") = res[1], 
        Rcpp::Named("Peri_C") = res[2], 
        Rcpp::Named("Peri_D") = res[3], 
        Rcpp::Named("Tele_C") = res[4], 
        Rcpp::Named("Tele_D") = res[5]
    );
}
