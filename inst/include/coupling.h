// Copyright (C) 2026 Wenbo Lyu
//
// This file is part of coupling.
//
// coupling is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// coupling is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with coupling. If not, see <https://www.gnu.org/licenses/>.

#ifndef COUPLING_COUPLING_H
#define COUPLING_COUPLING_H

// ============================================================
// Dependency Guard: Encourage best practices (non-blocking)
// ============================================================

#if defined(Rcpp_hpp) && !defined(COMPILING_COUPLING)
    #warning "It is recommended to include <coupling.h> alone, as it already includes <Rcpp.h>."
#endif

// ============================================================
// Core Dependencies (Auto-included for users)
// ============================================================

#include <Rcpp.h>

// ============================================================
// Module Headers (Organized by functionality)
// ============================================================

#include "coupling/numericutils.hpp"
#include "coupling/ccd.hpp"

// ============================================================
// Convenience Converters (Inline helpers for R/C++ interop)
// ============================================================

#include "coupling/convert.hpp"

#endif // COUPLING_COUPLING_H
