#ifndef COUPLING_METACOUPLING_HPP
#define COUPLING_METACOUPLING_HPP

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <RcppThread.h>
#include "coupling/ccd.hpp"

namespace coupling
{

namespace metacoupling
{

inline std::vector<double> metacoupling_c(
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::vector<double>>& swm_peri,
    const std::vector<std::vector<double>>& swm_tele,
    const std::string& method = "standard",
    size_t threads = 1
) {
    size_t n_units = mat.size();

    std::vector<double> result(3, std::numeric_limits<double>::quiet_NaN());
    if (n_units == 0) return result;
    
    if (threads <= 1) {
        for (size_t i = 0; i < n_units; ++i) {
            result[i] = ccd_c_single(mat[i], method);
        }
    } else {
        RcppThread::parallelFor(0, n_units, [&](size_t i) {
            result[i] = ccd_c_single(mat[i], method);
        }, threads);
    }

    return result;
}

inline std::vector<std::vector<double>> ccd(
    const std::vector<std::vector<double>>& mat,
    const std::vector<double>& weight,
    const std::string& method = "standard",
    size_t threads = 1
) {
    size_t n_units = mat.size(); // number of unit
    size_t p = mat[0].size(); // number of U values per unit

    std::vector<std::vector<double>> result(2, std::vector<double>(n_units, 0.0));
    
    std::vector<double> C_vals = ccd_c(mat, method, threads); // C values
    result[0] = C_vals;

    for (size_t i = 0; i < n_units; ++i) {
        double T_val = 0.0;
        for (size_t j = 0; j < p; ++j) {
            T_val += weight[j] * mat[i][j]; // T = \sum_{j=1}^n w_j \times U_j
        }
        result[1][i] = std::sqrt(C_vals[i] * T_val);
    }

    return result;
}

} // namespace ccd

}

#endif // COUPLING_METACOUPLING_HPP
