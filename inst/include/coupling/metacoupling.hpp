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
    size_t n_units = mat.size(); // number of observationa
    size_t p = mat[0].size(); // number of U values

    std::vector<std::vector<double>> result(
        n_units,
        std::vector<double>(3, std::numeric_limits<double>::quiet_NaN()));
    
    if (threads <= 1) {
        for (size_t s = 0; s < n_units; ++s) {
            result[s][0] = coupling::ccd::ccd_c_single(mat[s], method);
        }
    } else {
        RcppThread::parallelFor(0, n_units, [&](size_t s) {
            result[s][0] = coupling::ccd::ccd_c_single(mat[s], method);
        }, threads);
    }

    return result;
}

} // namespace metacoupling

}

#endif // COUPLING_METACOUPLING_HPP
