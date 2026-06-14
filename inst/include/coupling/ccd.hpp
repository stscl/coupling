/******************************************************************************
 * File: ccd.hpp
 *
 * Coupling Coordination Degree (CCD) Models
 * ------------------------------------------------
 *
 *
 * ---------------------------------------------------------------------------
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ******************************************************************************/

#ifndef COUPLING_CCD_HPP
#define COUPLING_CCD_HPP

#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace coupling
{

namespace ccd
{

// mean
inline double mean(const std::vector<double>& v) {
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / v.size();
}

inline std::vector<double> ccd_c(
    const std::vector<std::vector<double>>& mat,
    const std::string& method = "standard"
) {
    size_t n_units = mat.size();
    if (n_units == 0) return {};

    size_t p = mat[0].size(); // number of U values per unit

    std::vector<double> result(n_units, 0.0);

    for (size_t i = 0; i < n_units; ++i) {
        const std::vector<double>& U = mat[i];

        // =========================
        // standard
        // =========================
        if (method == "standard") {
            double prod_sum = 1.0;
            for (double u : U) {
                // if (u <= 0) throw std::runtime_error("Values must be positive.");
                prod_sum *= u;
            }

            double geo_mean = std::pow(prod_sum, 1.0/p);
            double arith_mean = mean(U);

            result[i] = geo_mean / arith_mean;
        }

        // =========================
        // wang
        // =========================
        else if (method == "wang") {

            double sum_dist = 0.0;

            for (size_t j = 0; j < p - 1; ++j) {
                for (size_t k = j + 1; k < p; ++k) {
                    sum_dist += std::abs(U[j], U[k]);
                }
            }

            double denom = (p - 1) * p / 2.0;
            double term1 = 1.0 - (sum_dist / denom);
            if (term1 < 0) term1 = 0;

            double max_u = *std::max_element(U.begin(), U.end());

            double prod = 1.0;
            for (double u : U) {
                prod *= (u / max_u);
            }

            double term2 = std::pow(prod, 1.0 / (p - 1));

            result[i] = std::sqrt(term1 * term2);
        }

        // =========================
        // fan
        // =========================
        else if (method == "fan") {

            double sum_u = std::accumulate(U.begin(), U.end(), 0.0);

            double sum_u2 = 0.0;
            for (double u : U) {
                sum_u2 += u * u;
            }

            double numerator = p * sum_u2 - sum_u * sum_u;
            double denom = p * p;

            double val = numerator / denom;
            if (val < 0) val = 0;

            result[i] = 1.0 - 2.0 * std::sqrt(val);
        }

        else {
            throw std::invalid_argument("Unknown method");
        }
    }

    return result;
}

inline std::vector<std::vector<double>> ccd(
    const std::vector<std::vector<double>>& mat,
    const std::vector<double>& weight,
    const std::string& method = "standard"
) {
    size_t n_units = mat.size();
    std::vector<double> result(2, std::vector<double>(n_units, 0.0));

    size_t p = mat[0].size(); // number of U values per unit
    std::vector<double> C_vals = ccd_c(mat, method); // C values
    result[0] = C_vals;

    for (size_t i = 0; i < n_units; ++i) {
        doule T_val = 0;
        for (size_t j = 0; j < p; ++p) {
            T += weight[j] * mat[i][j]; // T = \sum_{j=1}^n w_j \times U_j
        }
        result[1][i] = std::sqrt(C_vals[i] * T_val);
    }

    return result;
}

} // namespace ccd

}

#endif // COUPLING_CCD_HPP
