/******************************************************************************
 * File: metacoupling.hpp
 *
 * Metacoupling Model (Extension of CCD)
 * ------------------------------------------------
 *
 * This module implements the metacoupling framework, an extension of the
 * Coupling Coordination Degree (CCD) model, designed to measure both
 * intra-system and inter-system coupling relationships across spatial units.
 *
 * The metacoupling model integrates:
 *
 *      • intra-unit coupling (within a spatial unit)
 *      • peri-coupling (short-distance interactions)
 *      • tele-coupling (long-distance interactions)
 *
 * It further extends the CCD framework by introducing a permutation-based
 * cross-system coupling mechanism.
 *
 * ---------------------------------------------------------------------------
 * Mathematical formulation
 * ---------------------------------------------------------------------------
 *
 * (See https://github.com/stscl/coupling/discussions/8 for a better reading experience.)
 *
 * Let there be:
 *
 *      • n spatial units
 *      • p indicators per unit
 *
 * For unit i:
 *
 *      U_i = {U_{i1}, U_{i2}, ..., U_{ip}}
 *
 * ---------------------------------------------------------------------------
 * 1. Intra-system coupling (C_intra)
 * ---------------------------------------------------------------------------
 *
 * Standard CCD coupling computed for each unit:
 *
 *      C_intra(i) = CCD(U_i)
 *
 * where CCD is defined in ccd.hpp.
 *
 * ---------------------------------------------------------------------------
 * 2. Inter-system metacoupling (C_ij)
 * ---------------------------------------------------------------------------
 *
 * For two units i and j, metacoupling is defined by constructing
 * mixed systems through all possible indicator-level combinations:
 *
 *      V(mask) = {v_1, v_2, ..., v_p}
 *
 * where:
 *
 *      v_k = U_{ik}   if bit k in mask = 1
 *      v_k = U_{jk}   if bit k in mask = 0
 *
 * All binary masks are enumerated:
 *
 *      mask ∈ {0, 1}^p
 *
 * EXCLUDING:
 *
 *      • mask = 0            (all indicators from j)
 *      • mask = 2^p − 1      (all indicators from i)
 *
 * Total combinations used:
 *
 *      2^p − 2
 *
 * The inter-system coupling is:
 *
 *      C_ij = mean over all valid masks of CCD(V(mask))
 *
 * ---------------------------------------------------------------------------
 * 3. Spatial aggregation (peri / tele)
 * ---------------------------------------------------------------------------
 *
 * Given spatial weight matrices:
 *
 *      W_peri  (short-distance interactions)
 *      W_tele  (long-distance interactions)
 *
 * Both are assumed to be symmetric.
 *
 * Aggregated coupling:
 *
 *      C_peri(i) = Σ_j W_peri(i,j) * C_ij
 *      C_tele(i) = Σ_j W_tele(i,j) * C_ij
 *
 * ---------------------------------------------------------------------------
 * 4. Composite development index (T)
 * ---------------------------------------------------------------------------
 *
 * For each system:
 *
 *      T = Σ_k w_k * U_k
 *
 * For inter-system:
 *
 *      T_ij = mean over all permutations of T(V(mask))
 *
 * ---------------------------------------------------------------------------
 * 5. Coupling coordination degree (D)
 * ---------------------------------------------------------------------------
 *
 * Final metric:
 *
 *      D = sqrt(C × T)
 *
 * Applied to:
 *
 *      • intra-system
 *      • peri-coupling
 *      • tele-coupling
 *
 * ---------------------------------------------------------------------------
 * Input data format
 * ---------------------------------------------------------------------------
 *
 * mat : vector<vector<double>>
 *
 *      • rows = spatial units
 *      • columns = indicators
 *      • values MUST be normalized to [0,1]
 *
 * swm_peri / swm_tele:
 *
 *      • symmetric spatial weight matrices (n × n)
 *
 * weight:
 *
 *      • length = number of indicators
 *
  * ---------------------------------------------------------------------------
 * Output
 * ---------------------------------------------------------------------------
 *
 * metacoupling_c:
 *
 *      result is a 3 × n_units matrix:
 *
 *          result[0][i] = C_intra(i)
 *          result[1][i] = C_peri(i)
 *          result[2][i] = C_tele(i)
 *
 * metacoupling:
 *
 *      result is a 6 × n_units matrix:
 *
 *          result[0][i] = C_intra(i)
 *          result[1][i] = D_intra(i)
 *          result[2][i] = C_peri(i)
 *          result[3][i] = D_peri(i)
 *          result[4][i] = C_tele(i)
 *          result[5][i] = D_tele(i)
 *
 * ---------------------------------------------------------------------------
 * Notes
 * ---------------------------------------------------------------------------
 *
 * • The number of inter-system combinations grows exponentially (2^p − 2)
 *
 * • This framework captures cross-unit structural interactions
 *   at the indicator level, rather than treating systems as indivisible units
 *
 * • Input normalization is critical for meaningful interpretation
 *
 * • Computational cost increases rapidly with p
 *
 * ---------------------------------------------------------------------------
 * Author: Wenbo Lyu (Github: @SpatLyu)
 * License: GPL-3
 ******************************************************************************/

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
#include "coupling/numericutils.hpp"

namespace coupling
{

namespace metacoupling
{

inline std::vector<std::vector<double>> metacoupling_c( 
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::vector<double>>& swm_peri,
    const std::vector<std::vector<double>>& swm_tele,
    const std::string& method = "standard",
    size_t threads = 1
) {
    // Number of spatial units (observations)
    size_t n_units = mat.size();

    // Number of indicators (system dimensions)
    size_t p = mat[0].size();

    // Output:
    // [0][i] = intra-coupling (original CCD for unit i)
    // [1][i] = peri-coupling
    // [2][i] = tele-coupling
    std::vector<std::vector<double>> result(
        3,
        std::vector<double>(n_units, std::numeric_limits<double>::quiet_NaN())
    );

    if (p <= 1) return result;

    // Total number of binary combinations = 2^p
    size_t full_perm = 1ULL << p;

    // ============================================================
    // Compute metacoupling between unit i and j
    //
    // Definition:
    // - Enumerate ALL mixed systems
    // - Each dimension independently comes from i or j
    // - EXCLUDE:
    //      mask == 0           → all from j
    //      mask == full_perm-1 → all from i
    //
    // Total combinations used = 2^p - 2
    // ============================================================
    auto compute_pair = [&](size_t i, size_t j) {

        double sum_c = 0.0;
        size_t count = 0;

        for (size_t mask = 0; mask < full_perm; ++mask) {

            // ----------------------------------------------------
            // Remove trivial (non-mixed) systems
            // ----------------------------------------------------
            if (mask == 0 || mask == (full_perm - 1)) continue;

            // ----------------------------------------------------
            // Construct mixed system vector
            //
            // Bit k of mask determines source:
            //   1 → take from unit i
            //   0 → take from unit j
            //
            // Example (p = 3):
            // mask = 101 → [i, j, i]
            // ----------------------------------------------------
            std::vector<double> vec(p);

            for (size_t k = 0; k < p; ++k) {
                if (mask & (1ULL << k)) {
                    vec[k] = mat[i][k];
                } else {
                    vec[k] = mat[j][k];
                }
            }

            // Compute CCD coupling
            sum_c += coupling::ccd::ccd_c_single(vec, method);
            ++count;
        }

        // count should be exactly 2^p - 2
        return sum_c / static_cast<double>(count);
    };

    // ============================================================
    // Main worker for each spatial unit s
    // ============================================================
    auto worker = [&](size_t s) {

        // 1. Intra-coupling (original system)
        result[0][s] = coupling::ccd::ccd_c_single(mat[s], method);

        double peri_sum = 0.0;
        double tele_sum = 0.0;

        // 2. Iterate over all other units
        for (size_t j = 0; j < n_units; ++j) {

            double w_peri = swm_peri[s][j];
            double w_tele = swm_tele[s][j];

            // Skip if no interaction
            if (coupling::numericutils::doubleNearlyEqual(w_peri, 0.0) && 
                coupling::numericutils::doubleNearlyEqual(w_tele, 0.0)) continue;

            // Compute symmetric coupling between (s, j)
            double c_val = compute_pair(s, j);

            // 3. Weighted aggregation
            // Since SWM is symmetric, no need to compute (j, s) separately
            if (!coupling::numericutils::doubleNearlyEqual(w_peri, 0.0)) {
                peri_sum += w_peri * c_val;
            }

            if (!coupling::numericutils::doubleNearlyEqual(w_tele, 0.0)) {
                tele_sum += w_tele * c_val;
            }
        }

        result[1][s] = peri_sum;
        result[2][s] = tele_sum;
    };
    
    // ============================================================
    // Check for available threads
    // ============================================================
    if (threads == 0) threads = 1;
    size_t hw = std::thread::hardware_concurrency();
    if (hw > 0) threads = std::min(threads, hw);

    // ============================================================
    // Parallel or serial execution
    // ============================================================
    if (threads <= 1) {
        for (size_t s = 0; s < n_units; ++s) {
            worker(s);
        }
    } else {
        RcppThread::parallelFor(0, n_units, worker, threads);
    }

    return result;
}

inline std::vector<std::vector<double>> metacoupling( 
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::vector<double>>& swm_peri,
    const std::vector<std::vector<double>>& swm_tele,
    const std::vector<double>& weight,
    const std::string& method = "standard",
    size_t threads = 1
) {
    size_t n_units = mat.size();
    size_t p = mat[0].size();
    
    // ============================================================
    // Output: 6 × n_units
    //
    // [0] = Intra_C
    // [1] = Intra_D
    // [2] = Peri_C
    // [3] = Peri_D
    // [4] = Tele_C
    // [5] = Tele_D
    // ============================================================
    std::vector<std::vector<double>> result(
        6,
        std::vector<double>(n_units, std::numeric_limits<double>::quiet_NaN())
    );

    if (p <= 1) return result;

    size_t full_perm = 1ULL << p;

    // ============================================================
    // Check for available threads
    // ============================================================
    if (threads == 0) threads = 1;
    size_t hw = std::thread::hardware_concurrency();
    if (hw > 0) threads = std::min(threads, hw);

    // ============================================================
    // Step 1: reuse existing C computation
    // ============================================================
    auto C_res = metacoupling_c(mat, swm_peri, swm_tele, method, threads);

    // ============================================================
    // Step 2: T computation
    // ============================================================
    auto compute_T_intra = [&](size_t i) {
        double t = 0.0;
        for (size_t k = 0; k < p; ++k) {
            t += weight[k] * mat[i][k];
        }
        return t;
    };

    auto compute_pair_T = [&](size_t i, size_t j) {

        double sum_t = 0.0;
        size_t count = 0;

        for (size_t mask = 0; mask < full_perm; ++mask) {

            if (mask == 0 || mask == (full_perm - 1)) continue;

            double t_val = 0.0;

            for (size_t k = 0; k < p; ++k) {
                if (mask & (1ULL << k)) {
                    t_val += weight[k] * mat[i][k];
                } else {
                    t_val += weight[k] * mat[j][k];
                }
            }

            sum_t += t_val;
            ++count;
        }

        return sum_t / static_cast<double>(count);
    };

    auto worker4D = [&](size_t s) {
        // ---- intra ----
        double C_intra = C_res[0][s];
        double T_intra = compute_T_intra(s);

        result[0][s] = C_intra;
        result[1][s] = std::sqrt(C_intra * T_intra);

        double peri_T_sum = 0.0;
        double tele_T_sum = 0.0;

        // ---- interactions ----
        for (size_t j = 0; j < n_units; ++j) {

            double w_peri = swm_peri[s][j];
            double w_tele = swm_tele[s][j];

            if (coupling::numericutils::doubleNearlyEqual(w_peri, 0.0) && 
                coupling::numericutils::doubleNearlyEqual(w_tele, 0.0)) continue;

            double t_val = compute_pair_T(s, j);

            if (!coupling::numericutils::doubleNearlyEqual(w_peri, 0.0)) {
                peri_T_sum += w_peri * t_val;
            }

            if (!coupling::numericutils::doubleNearlyEqual(w_tele, 0.0)) {
                tele_T_sum += w_tele * t_val;
            }
        }

        double C_peri = C_res[1][s];
        double C_tele = C_res[2][s];

        result[2][s] = C_peri;
        result[3][s] = std::sqrt(C_peri * peri_T_sum);

        result[4][s] = C_tele;
        result[5][s] = std::sqrt(C_tele * tele_T_sum);
    };

    // ============================================================
    // Step 3: assemble results
    // ============================================================
    if (threads <= 1) {
        for (size_t s = 0; s < n_units; ++s) {
            worker4D(s);
        }
    } else {
        RcppThread::parallelFor(0, n_units, worker4D, threads);
    }

    return result;
}

} // namespace metacoupling

}

#endif // COUPLING_METACOUPLING_HPP
