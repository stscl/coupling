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
    // [i][0] = intra-coupling (original CCD for unit i)
    // [i][1] = peri-coupling
    // [i][2] = tele-coupling
    std::vector<std::vector<double>> result(
        n_units,
        std::vector<double>(3, std::numeric_limits<double>::quiet_NaN())
    );

    // Total number of full binary combinations = 2^p
    // Each bit represents whether a dimension comes from unit i or unit j
    size_t full_perm = 1ULL << p;

    // ============================================================
    // Compute coupling between a pair (i, j)
    // using HALF of the full permutation space: 2^(p-1)
    // ============================================================
    auto compute_pair = [&](size_t i, size_t j) {

        double sum_c = 0.0;
        size_t count = 0;

        // Iterate over all 2^p binary masks
        for (size_t mask = 0; mask < full_perm; ++mask) {

            // ----------------------------------------------------
            // Symmetry reduction:
            // For any mask, its complement (~mask) produces
            // the same combination but with i and j swapped.
            //
            // Example:
            // mask      = 0101
            // complement= 1010
            //
            // These two correspond to equivalent mixed systems.
            //
            // To avoid double counting, we only keep one of them:
            // enforce: mask <= complement
            // ----------------------------------------------------
            size_t complement = (~mask) & (full_perm - 1);

            if (mask > complement) continue;

            // Construct a mixed indicator vector
            std::vector<double> vec(p);

            // ----------------------------------------------------
            // Bitwise construction:
            // For each dimension k:
            //
            // if bit k == 1 → take value from unit i
            // if bit k == 0 → take value from unit j
            //
            // This ensures:
            // - all possible combinations are explored
            // - no dimension is artificially fixed
            // ----------------------------------------------------
            for (size_t k = 0; k < p; ++k) {
                if (mask & (1ULL << k)) {
                    vec[k] = mat[i][k];
                } else {
                    vec[k] = mat[j][k];
                }
            }

            // Compute CCD coupling for this mixed system
            sum_c += coupling::ccd::ccd_c_single(vec, method);
            ++count;
        }

        // The number of retained combinations is exactly 2^(p-1)
        return sum_c / static_cast<double>(count);
    };

    // ============================================================
    // Main worker for each spatial unit s
    // ============================================================
    auto worker = [&](size_t s) {

        // 1. Intra-coupling (original system)
        result[s][0] = coupling::ccd::ccd_c_single(mat[s], method);

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

        result[s][1] = peri_sum;
        result[s][2] = tele_sum;
    };

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

} // namespace metacoupling

}

#endif // COUPLING_METACOUPLING_HPP
