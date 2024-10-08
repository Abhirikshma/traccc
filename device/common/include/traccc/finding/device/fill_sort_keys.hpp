/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_candidate.hpp"

namespace traccc::device {

/// Function used for fill key container
///
/// @param[in] globalIndex   The index of the current thread
/// @param[in] params_view   The input parameters
/// @param[out] keys_view    The key values
/// @param[out] ids_view     The param ids
///
TRACCC_HOST_DEVICE inline void fill_sort_keys(
    std::size_t globalIndex,
    bound_track_parameters_collection_types::const_view params_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/fill_sort_keys.ipp"
