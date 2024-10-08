/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/cell.hpp"
#include "traccc/edm/seed.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/edm/track_candidate.hpp"
#include "traccc/io/data_format.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"

// System include(s).
#include <cstddef>
#include <string_view>

namespace traccc::io {

/// Function for cell file writing
///
/// @param event is the event index
/// @param directory is the directory for the output cell file
/// @param format is the data format (e.g. csv or binary) of output file
/// @param cells is the cell collection to write
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           traccc::cell_collection_types::const_view cells);

/// Function for hit file writing
///
/// @param event is the event index
/// @param directory is the directory for the output spacepoint file
/// @param format is the data format (e.g. csv or binary) of output file
/// @param spacepoints is the spacepoint collection to write
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           spacepoint_collection_types::const_view spacepoints);

/// Function for measurement file writing
///
/// @param event is the event index
/// @param directory is the directory for the output measurement file
/// @param format is the data format (e.g. csv or binary) of output file
/// @param measurements is the measurement collection to write
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           measurement_collection_types::const_view measurements);

/// Function for seed writing
///
/// @param event is the event index
/// @param directory is the directory for the output seed file
/// @param format is the data format (obj only right now) of output file
/// @param seeds is the seed collection to write
/// @param spacepoints is the spacepoint collection the seeds are made of
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format, seed_collection_types::const_view seeds,
           spacepoint_collection_types::const_view spacepoints);

/// Function for track candidate writing
///
/// @param event is the event index
/// @param directory is the directory for the output seed file
/// @param format is the data format (obj only right now) of output file
/// @param tracks is the track candidate container to write
/// @param detector is the Detray detector describing the geometry
///
void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           track_candidate_container_types::const_view tracks,
           const detray::detector<>& detector);

}  // namespace traccc::io
