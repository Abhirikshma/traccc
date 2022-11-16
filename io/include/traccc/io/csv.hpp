/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/cell.hpp"
#include "traccc/edm/cluster.hpp"
#include "traccc/edm/container.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/geometry/geometry.hpp"
#include "traccc/io/detail/json_digitization_config.hpp"

// Acts include(s)
#include "Acts/Geometry/GeometryHierarchyMap.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// DFE include(s).
#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>

// System include(s).
#include <cassert>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

/// reader
namespace traccc {

/// reader

struct csv_cell {

    uint64_t geometry_id = 0;
    uint64_t hit_id = 0;
    channel_id channel0 = 0;
    channel_id channel1 = 0;
    scalar timestamp = 0.;
    scalar value = 0.;

    // geometry_id,hit_id,channel0,channel1,timestamp,value
    DFE_NAMEDTUPLE(csv_cell, geometry_id, hit_id, channel0, channel1, timestamp,
                   value);
};

using cell_reader = dfe::NamedTupleCsvReader<csv_cell>;

struct csv_meas_hit_id {

    uint64_t measurement_id = 0;
    uint64_t hit_id = 0;

    // geometry_id,hit_id,channel0,channel1,timestamp,value
    DFE_NAMEDTUPLE(csv_meas_hit_id, measurement_id, hit_id);
};

using meas_hit_id_reader = dfe::NamedTupleCsvReader<csv_meas_hit_id>;

struct csv_fatras_hit {
    uint64_t particle_id = 0;
    uint64_t geometry_id = 0;
    scalar tx = 0;
    scalar ty = 0;
    scalar tz = 0;
    scalar tt = 0;
    scalar tpx = 0;
    scalar tpy = 0;
    scalar tpz = 0;
    scalar te = 0;
    scalar deltapx = 0;
    scalar deltapy = 0;
    scalar deltapz = 0;
    scalar deltae = 0;
    uint64_t index = 0;

    DFE_NAMEDTUPLE(csv_fatras_hit, particle_id, geometry_id, tx, ty, tz, tt,
                   tpx, tpy, tpz, te, deltapx, deltapy, deltapz, deltae, index);
};

using fatras_hit_reader = dfe::NamedTupleCsvReader<csv_fatras_hit>;

struct csv_particle {
    uint64_t particle_id = 0;
    int particle_type = 0;
    int process = 0;
    scalar vx;
    scalar vy;
    scalar vz;
    scalar vt;
    scalar px;
    scalar py;
    scalar pz;
    scalar m;
    scalar q;

    DFE_NAMEDTUPLE(csv_particle, particle_id, particle_type, process, vx, vy,
                   vz, vt, px, py, pz, m, q);
};

using fatras_particle_reader = dfe::NamedTupleCsvReader<csv_particle>;

/// writer

struct csv_measurement {

    uint64_t geometry_id = 0;
    std::string local_key = "";
    scalar local0 = 0.;
    scalar local1 = 0.;
    scalar phi = 0.;
    scalar theta = 0.;
    scalar time = 0.;
    scalar var_local0 = 0.;
    scalar var_local1 = 0.;
    scalar var_phi = 0.;
    scalar var_theta = 0.;
    scalar var_time = 0.;

    DFE_NAMEDTUPLE(csv_measurement, geometry_id, local0, local1, phi, theta,
                   time, var_local0, var_local1, var_phi, var_theta, var_time);
};

using measurement_reader = dfe::NamedTupleCsvReader<csv_measurement>;
using measurement_writer = dfe::NamedTupleCsvWriter<csv_measurement>;

struct csv_spacepoint {

    uint64_t geometry_id = 0;
    scalar x, y, z;

    DFE_NAMEDTUPLE(csv_spacepoint, geometry_id, x, y, z);
};

using spacepoint_writer = dfe::NamedTupleCsvWriter<csv_spacepoint>;

struct csv_seed {
    scalar weight;
    scalar z_vertex;
    scalar x_b, y_b, z_b;
    scalar varR_b, varZ_b;
    scalar x_m, y_m, z_m;
    scalar varR_m, varZ_m;
    scalar x_t, y_t, z_t;
    scalar varR_t, varZ_t;

    DFE_NAMEDTUPLE(csv_seed, weight, z_vertex, x_b, y_b, z_b, varR_b, varZ_b,
                   x_m, y_m, z_m, varR_m, varZ_m, x_t, y_t, z_t, varR_t,
                   varZ_t);
};

using seed_writer = dfe::NamedTupleCsvWriter<csv_seed>;

struct csv_bound_track_parameters {
    scalar loc0;
    scalar loc1;
    scalar theta;
    scalar phi;
    scalar qoverp;
    scalar time;
    scalar cov00;
    scalar cov01, cov11;
    scalar cov02, cov12, cov22;
    scalar cov03, cov13, cov23, cov33;
    scalar cov04, cov14, cov24, cov34, cov44;
    scalar cov05, cov15, cov25, cov35, cov45, cov55;
    unsigned int surface_id;

    DFE_NAMEDTUPLE(csv_bound_track_parameters, loc0, loc1, theta, phi, qoverp,
                   time, cov00, cov01, cov11, cov02, cov12, cov22, cov03, cov13,
                   cov23, cov33, cov04, cov14, cov24, cov34, cov44, cov05,
                   cov15, cov25, cov35, cov45, cov55, surface_id);
};

using bound_track_parameters_writer =
    dfe::NamedTupleCsvWriter<csv_bound_track_parameters>;

struct cell_counter {
    geometry_id module = 0;
    unsigned int nCells = 0;
};  // struct cell_counter

/// Declare all cell counter collection types
using cell_counter_collection_types = collection_types<cell_counter>;

/// Read the collection of cells per module and fill into a collection
///
/// @param creader The cellreader type
/// @param resource The memory resource to use for the return value
/// @param tfmap the (optional) transform map
/// @param digi_cfg the (optional) digitization configuration
/// @param max_cells the (optional) maximum number of cells to be read in
inline cell_container_types::host read_cells(
    cell_reader& creader, vecmem::memory_resource& resource,
    const traccc::geometry* tfmap = nullptr,
    const Acts::GeometryHierarchyMap<digitization_config>* digi_cfg = nullptr,
    unsigned int max_cells = std::numeric_limits<unsigned int>::max()) {

    // Create cell counter vector
    cell_counter_collection_types::host cell_counts(&resource);
    // Create cell vector
    cell_collection_types::host allCells(&resource);
    allCells.reserve(50000);
    cell_counts.reserve(5000);

    unsigned int read_cells = 0;
    csv_cell iocell;

    // Read all cells from input file
    while (creader.read(iocell)) {

        cell c = cell{iocell.channel0, iocell.channel1, iocell.value,
                      iocell.timestamp};

        // Add current cell to vector
        allCells.push_back(c);

        geometry_id id = iocell.geometry_id;

        auto rit = std::find_if(
            cell_counts.rbegin(), cell_counts.rend(),
            [&id](const cell_counter& cc) { return cc.module == id; });

        if (rit == cell_counts.rend()) {
            // Add module to cell counter vector
            cell_counts.push_back({id, 1});
        } else {
            // Increment number of cells for this module
            ++(rit->nCells);
        }

        if (++read_cells >= max_cells) {
            break;
        }
    }

    const std::size_t size = cell_counts.size();

    // Construct the result container and reserve sizes for jagged vector.
    cell_container_types::host result(size, &resource);
    for (std::size_t i = 0; i < size; ++i) {
        result.get_items().at(i).reserve(cell_counts[i].nCells);
    }

    std::size_t idx = 0;
    // Fill result container
    for (std::size_t i = 0; i < size; ++i) {
        cell_module module;
        module.module = cell_counts[i].module;

        if (tfmap != nullptr) {
            if (!tfmap->contains(module.module)) {
                throw std::runtime_error(
                    "Could not find placement for geometry ID " +
                    std::to_string(module.module));
            }
            module.placement = (*tfmap)[module.module];
        }

        if (digi_cfg != nullptr) {
            Acts::GeometryIdentifier id(module.module);

            auto geo_it = digi_cfg->find(id);
            if (geo_it == digi_cfg->end()) {
                throw std::runtime_error(
                    "Could not find digitization config for geometry ID " +
                    std::to_string(module.module));
            }

            const auto& binning_data = geo_it->segmentation.binningData();

            module.pixel =
                pixel_data{binning_data[0].min, binning_data[1].min,
                           binning_data[0].step, binning_data[1].step};
        }

        for (unsigned int j = 0; j < cell_counts[i].nCells; ++j) {
            const cell& c = allCells[idx++];
            // Fill item for this cell
            result.get_items().at(i).push_back(c);
            module.range0[0] = std::min(module.range0[0], c.channel0);
            module.range0[1] = std::max(module.range0[1], c.channel0);
            module.range1[0] = std::min(module.range1[0], c.channel1);
            module.range1[1] = std::max(module.range1[1], c.channel1);
        }
        // Fill header for this module
        result.get_headers().at(i) = module;

        auto& cells = result.get_items().at(i);
        std::sort(cells.begin(), cells.end(), [](const auto& a, const auto& b) {
            return a.channel1 < b.channel1;
        });
    }
    return result;
}

/// Read the collection of hits per module and fill into a collection
///
/// @param hreader The hit reader type
/// @param resource The memory resource to use for the return value
inline spacepoint_container_types::host read_hits(
    fatras_hit_reader& hreader, vecmem::memory_resource& resource,
    const traccc::geometry* tfmap = nullptr,
    unsigned int max_hits = std::numeric_limits<unsigned int>::max()) {
    spacepoint_container_types::host result(&resource);

    unsigned int read_hits = 0;
    csv_fatras_hit iohit;

    while (hreader.read(iohit)) {
        geometry_id geom_id = iohit.geometry_id;
        auto placement = (*tfmap)[geom_id];

        point3 position({iohit.tx, iohit.ty, iohit.tz});
        auto local = placement.point_to_local(position);
        measurement m({point2({local[0], local[1]}), variance2({0., 0.})});
        spacepoint sp({position, m});

        const spacepoint_container_types::host::header_vector& headers =
            result.get_headers();

        auto rit = std::find(headers.rbegin(), headers.rend(), geom_id);

        if (rit == headers.rend()) {
            result.push_back(geom_id, vecmem::vector<spacepoint>({sp}));
        } else {
            // The reverse iterator.base() returns the equivalent normal
            // iterator shifted by 1, so that the (r)end and (r)begin iterators
            // match consistently, due to the extra past-the-last element
            auto idx = std::distance(headers.begin(), rit.base()) - 1;
            result.at(idx).items.push_back(sp);
        }

        if (++read_hits >= max_hits) {
            break;
        }
    }

    return result;
}

/// Read the collection of measurements per module and fill into a collection
///
/// @param mreader The measurement reader type
/// @param resource The memory resource to use for the return value
inline measurement_container_types::host read_measurements(
    measurement_reader& mreader, vecmem::memory_resource& resource,
    unsigned int max_measurements = std::numeric_limits<unsigned int>::max()) {
    measurement_container_types::host result(&resource);

    unsigned int read_measurements = 0;
    csv_measurement io_measurement;

    while (mreader.read(io_measurement)) {
        cell_module module;

        module.module = io_measurement.geometry_id;

        measurement meas;
        meas.local = {io_measurement.local0, io_measurement.local1};
        meas.variance = {io_measurement.var_local0, io_measurement.var_local1};

        const measurement_container_types::host::header_vector& headers =
            result.get_headers();

        auto rit = std::find(headers.rbegin(), headers.rend(), module);

        if (rit == headers.rend()) {
            result.push_back(module, vecmem::vector<measurement>({meas}));
        } else {
            // The reverse iterator.base() returns the equivalent normal
            // iterator shifted by 1, so that the (r)end and (r)begin iterators
            // match consistently, due to the extra past-the-last element
            auto idx = std::distance(headers.begin(), rit.base()) - 1;
            result.at(idx).items.push_back(meas);
        }

        if (++read_measurements >= max_measurements) {
            break;
        }
    }
    return result;
}

}  // namespace traccc
