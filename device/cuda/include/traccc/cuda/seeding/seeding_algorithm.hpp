/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/cuda/seeding/details/seed_finding.hpp"
#include "traccc/cuda/seeding/details/spacepoint_binning.hpp"
#include "traccc/cuda/utils/stream.hpp"

// Project include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <memory>

namespace traccc::cuda {

/// Main algorithm for performing the track seeding on an NVIDIA GPU
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class seeding_algorithm : public algorithm<edm::seed_collection::buffer(
                              const edm::spacepoint_collection::const_view&)>,
                          public messaging {

    public:
    /// Constructor for the seed finding algorithm
    ///
    /// @param mr The memory resource to use
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param str The CUDA stream to perform the operations in
    ///
    seeding_algorithm(
        const seedfinder_config& finder_config,
        const spacepoint_grid_config& grid_config,
        const seedfilter_config& filter_config,
        const traccc::memory_resource& mr, vecmem::copy& copy, stream& str,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Operator executing the algorithm.
    ///
    /// @param spacepoints is a view of all spacepoints in the event
    /// @return the buffer of track seeds reconstructed from the spacepoints
    ///
    output_type operator()(const edm::spacepoint_collection::const_view&
                               spacepoints) const override;

    private:
    /// Tool performing the spacepoint binning
    details::spacepoint_binning m_binning;
    /// Tool performing the seed finding
    details::seed_finding m_finding;

};  // class seeding_algorithm

}  // namespace traccc::cuda
