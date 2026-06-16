#ifndef COUPLING_METACOUPLING_HPP
#define COUPLING_METACOUPLING_HPP

#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <RcppThread.h>
#include "coupling/ccd.hpp"

namespace coupling
{

namespace metacoupling
{

// mean
inline double mean(const std::vector<double>& v) {
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    return s / v.size();
}

inline double ccd_c_single(
    const std::vector<double>& vec,
    const std::string& method = "standard"
) {
    size_t p = vec.size(); // number of U values
    
    double C_val = 0.0;

    // =========================
    // standard
    // =========================
    if (method == "standard") {
        double prod_sum = 1.0;
        for (double u : vec) {
            // if (u <= 0) throw std::runtime_error("Values must be positive.");
            prod_sum *= u;
        }

        double geo_mean = std::pow(prod_sum, 1.0/p);
        double arith_mean = mean(vec);
        
        if (coupling::numericutils::doubleNearlyEqual(arith_mean, 0.0)) {
            return 0.0;
        }

        C_val = geo_mean / arith_mean;
    }

    // =========================
    // wang
    // =========================
    else if (method == "wang") {
        double sum_dist = 0.0;

        for (size_t i = 0; i < p - 1; ++i) {
            for (size_t j = i + 1; j < p; ++j) {
                sum_dist += std::abs(vec[i] - vec[j]);
            }
        }

        double denom = (p - 1) * p / 2.0;
        double term1 = 1.0 - (sum_dist / denom);
        // if (term1 < 0) term1 = 0;

        double max_u = *std::max_element(vec.begin(), vec.end());

        if (coupling::numericutils::doubleNearlyEqual(max_u, 0.0)) {
            return 0.0;
        }

        double prod = 1.0;
        for (double u : vec) {
            prod *= (u / max_u);
        }

        double term2 = std::pow(prod, 1.0 / (p - 1));

        C_val = std::sqrt(term1 * term2);
    }

    // =========================
    // fan
    // =========================
    else if (method == "fan") {
        double sum_u = std::accumulate(vec.begin(), vec.end(), 0.0);
        // if (coupling::numericutils::doubleNearlyEqual(sum_u, 0.0)) {
        //     return 0.0;
        // }

        double sum_u2 = 0.0;
        for (double u : vec) {
            sum_u2 += u * u;
        }
        // if (coupling::numericutils::doubleNearlyEqual(sum_u2, 0.0)) {
        //     return 0.0;
        // }

        double numerator = p * sum_u2 - sum_u * sum_u;
        double denom = p * p;

        double val = numerator / denom;
        if (val < 0) val = 0;

        C_val = 1.0 - 2.0 * std::sqrt(val);
    }

    else {
        throw std::invalid_argument("Unknown method");
    }
    
    return std::clamp(C_val, 0.0, 1.0);
}

inline std::vector<double> ccd_c(
    const std::vector<std::vector<double>>& mat,
    const std::string& method = "standard",
    size_t threads = 1
) {
    size_t n_units = mat.size();
    if (n_units == 0) return {};

    std::vector<double> result(n_units, 0.0);
    
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
