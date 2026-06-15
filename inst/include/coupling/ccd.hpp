/******************************************************************************
 * File: ccd.hpp
 *
 * Coupling Coordination Degree (CCD) Models
 * ------------------------------------------------
 *
 * This module implements the Coupling Coordination Degree (CCD) framework,
 * a widely used approach for measuring the interaction, coupling strength,
 * and coordinated development among multiple subsystems or indicators.
 *
 * The CCD model evaluates the internal consistency of a system composed of
 * multiple variables and further quantifies its coordinated development level
 * by combining coupling strength with a composite development index.
 *
 * ---------------------------------------------------------------------------
 * Mathematical formulation
 * ---------------------------------------------------------------------------
 * 
 * (See https://github.com/stscl/coupling/discussions/3 for a better reading experience.)
 *
 * Let a system consist of n indicators:
 *
 *      U = {U1, U2, ..., Un}
 *
 * For each spatial unit (or observation), the CCD model is computed based on
 * the vector U.
 *
 * ---------------------------------------------------------------------------
 * 1. Coupling degree (C)
 * ---------------------------------------------------------------------------
 *
 * Three alternative formulations are supported:
 *
 * (1) Standard coupling model:
 *
 *      C = [ (∏_{i=1}^n U_i) / ( (1/n * ∑_{i=1}^n U_i)^n ) ]^(1/n)
 *
 * This formulation reflects the balance among indicators using the ratio of
 * geometric mean to arithmetic mean.
 *
 * ---------------------------------------------------------------------------
 *
 * (2) Wang et al. model:
 *
 *      C = sqrt(
 *              [ 1 − ( ∑_{i>j, j = 1,2,...n-1} |U_i − U_j| ) / ( ∑_{m=1}^{n-1} m ) ]
 *              × ( ∏_{i=1}^n (U_i / max(U)) )^(1/(n−1))
 *          )
 *
 * This formulation incorporates both:
 *
 *      • dispersion among indicators (pairwise differences)
 *      • relative development level (normalized by max value)
 *
 * ---------------------------------------------------------------------------
 *
 * (3) Fan et al. model:
 *
 *      C = 1 − 2 * sqrt(
 *              ( n * ∑_{i=1}^n U_i^2 − (∑_{i=1}^n U_i)^2 ) / n^2
 *          )
 *
 * This formulation is based on variance structure and reflects inequality
 * among indicators.
 *
 * ---------------------------------------------------------------------------
 * 2. Composite development index (T)
 * ---------------------------------------------------------------------------
 *
 * A weighted linear aggregation of indicators:
 *
 *      T = ∑_{i=1}^n w_i * U_i
 *
 * subject to:
 *
 *      ∑_{i=1}^n w_i = 1
 *
 * where w_i are user-provided weights.
 *
 * ---------------------------------------------------------------------------
 * 3. Coupling coordination degree (D)
 * ---------------------------------------------------------------------------
 *
 * The final CCD index is defined as:
 *
 *      D = sqrt(C × T)
 *
 * This combines system interaction (C) with development level (T),
 * providing an overall measure of coordinated development.
 *
 * ---------------------------------------------------------------------------
 * Input data format
 * ---------------------------------------------------------------------------
 *
 * The input data is organized as:
 *
 *      mat : matrix (vector<vector<double>>)
 *
 * where:
 *
 *      • Each row corresponds to one spatial unit (or observation)
 *      • Each column corresponds to one indicator U_i
 *
 * That is:
 *
 *      mat[i] = {U_1, U_2, ..., U_n} for spatial unit i
 *
 * Weights:
 *
 *      weight : vector<double>, length = number of indicators
 *
 * ---------------------------------------------------------------------------
 * Output
 * ---------------------------------------------------------------------------
 *
 * The main functions return:
 *
 *      • C values: vector<double> by ccd_c
 *          Coupling degree for each spatial unit
 *
 *      • CCD values: vector<vector<double>> by ccd
 *          result[0] = C values
 *          result[1] = D values
 *
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * • All computations are performed independently for each spatial unit.
 *
 * • Input values are assumed to be non-negative. For the standard model,
 *   strictly positive values are required due to the geometric mean.
 *
 * • Numerical safeguards are applied to avoid invalid operations such as
 *   negative values under square roots caused by floating-point errors.
 *
 * • The Wang formulation uses pairwise absolute differences as a measure
 *   of dispersion.
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
                    sum_dist += std::abs(U[j] - U[k]);
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
    size_t n_units = mat.size(); // number of unit
    size_t p = mat[0].size(); // number of U values per unit

    std::vector<std::vector<double>> result(2, std::vector<double>(n_units, 0.0));
    
    std::vector<double> C_vals = ccd_c(mat, method); // C values
    result[0] = C_vals;

    for (size_t i = 0; i < n_units; ++i) {
        doule T_val = 0.0;
        for (size_t j = 0; j < p; ++j) {
            T_val += weight[j] * mat[i][j]; // T = \sum_{j=1}^n w_j \times U_j
        }
        result[1][i] = std::sqrt(C_vals[i] * T_val);
    }

    return result;
}

} // namespace ccd

}

#endif // COUPLING_CCD_HPP
