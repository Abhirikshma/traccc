/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/cell.hpp"
#include "traccc/geometry/geometry.hpp"
#include "traccc/io/digitization_config.hpp"
#include "traccc/io/reader_edm.hpp"

// System include(s).
#include <cstdint>
#include <map>
#include <string_view>

namespace traccc::io::csv {

/// Read cell information from a specific CSV file
///
/// @param out A cell (host) collection & a cell_module collection
/// @param filename The file to read the cell data from
/// @param geom The description of the detector geometry
/// @param dconfig The detector's digitization configuration
/// @param bardoce_map An object to perform barcode re-mapping with
/// @param deduplicate Whether to deduplicate the cells
///
void read_cells(cell_reader_output& out, std::string_view filename,
                const geometry* geom = nullptr,
                const digitization_config* dconfig = nullptr,
                const std::map<std::uint64_t, detray::geometry::barcode>*
                    barcode_map = nullptr,
                bool deduplicate = true);

}  // namespace traccc::io::csv
