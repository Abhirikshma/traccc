/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_candidate.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/fitting/triplet_fit/triplet_fitter.hpp"
#include "traccc/utils/algorithm.hpp"

#include <iostream>

namespace traccc {

/// Fitting algorithm for a set of tracks
template <typename fitter_t>
class fitting_algorithm
    : public algorithm<track_state_container_types::host(
          const typename fitter_t::detector_type&,
          const typename fitter_t::bfield_type&,
          const typename track_candidate_container_types::host&)> {

    public:
    using algebra_type = typename fitter_t::algebra_type;
    using bfield_type = typename fitter_t::bfield_type;
    /// Configuration type
    using config_type = typename fitter_t::config_type;

    /// Constructor for the fitting algorithm
    ///
    /// @param cfg  Configuration object
    fitting_algorithm(const config_type& cfg) : m_cfg(cfg) {}

    /// Run the algorithm
    ///
    /// @param track_candidates the candidate measurements from track finding
    /// @return the container of the fitted track parameters
    track_state_container_types::host operator()(
        const typename fitter_t::detector_type& det,
        const typename fitter_t::bfield_type& field,
        const typename track_candidate_container_types::host& track_candidates)
        const override {

        fitter_t fitter(det, field, m_cfg);

        track_state_container_types::host output_states;

        // The number of tracks
        std::size_t n_tracks = track_candidates.size();

        // Iterate over tracks
        for (std::size_t i = 0; i < n_tracks; i++) {

            // Seed parameter
            const auto& seed_param = track_candidates[i].header;

            // Make a vector of track state
            auto& cands = track_candidates[i].items;
            vecmem::vector<track_state<algebra_type>> input_states;
            input_states.reserve(cands.size());
            for (auto& cand : cands) {
                input_states.emplace_back(cand);
            }

            // Make a fitter state
            typename fitter_t::state fitter_state(input_states);

            // Run fitter
            fitter.fit(seed_param, fitter_state);

            output_states.push_back(
                std::move(fitter_state.m_fit_res),
                std::move(fitter_state.m_fit_actor_state.m_track_states));
        }

        return output_states;
    }

    /// Config object
    config_type m_cfg;
};

/// Partial specialization for Triplet Fitter
template<typename stepper_t, typename navigator_t>
class fitting_algorithm<traccc::triplet_fitter<stepper_t, navigator_t>>
    : public algorithm<track_state_container_types::host(
          const typename traccc::triplet_fitter<stepper_t, navigator_t>::detector_type&,
          const typename traccc::triplet_fitter<stepper_t, navigator_t>::bfield_type&,
          const typename track_candidate_container_types::host&)> {

    public:
    using algebra_type = typename traccc::triplet_fitter<stepper_t, navigator_t>::algebra_type;
    using bfield_type = typename traccc::triplet_fitter<stepper_t, navigator_t>::bfield_type;
    /// Configuration type
    using config_type = typename traccc::triplet_fitter<stepper_t, navigator_t>::config_type;

    /// Constructor for the fitting algorithm
    ///
    /// @param cfg  Configuration object
    fitting_algorithm(const config_type& cfg) : m_cfg(cfg) {}

    /// Run the algorithm
    ///
    /// @param track_candidates the candidate measurements from track finding
    /// @return the container of the fitted track parameters
    track_state_container_types::host operator()(
        const typename traccc::triplet_fitter<stepper_t, navigator_t>::detector_type& det,
        const typename traccc::triplet_fitter<stepper_t, navigator_t>::bfield_type& field,
        const typename track_candidate_container_types::host& track_candidates)
        const override {

        traccc::triplet_fitter<stepper_t, navigator_t> fitter(det, field, m_cfg);

        track_state_container_types::host output_states;

        // The number of tracks
        std::size_t n_tracks = track_candidates.size();

        // Iterate over tracks
        for (std::size_t i = 0; i < n_tracks; i++) {

            std::cout << "\nFitting track # " << i << std::endl;

            // Seed parameter
            // const auto& seed_param = track_candidates[i].header;

            // Make a vector of track state
            auto& cands = track_candidates[i].items;

            // std::cout << cands.size() << " measurements in this track" << std::endl;

            vecmem::vector<track_state<algebra_type>> input_states;
            input_states.reserve(cands.size());
            for (auto& cand : cands) {
                input_states.emplace_back(cand);
                // std::cout << "\n\tlocal position of meas: " << cand.local[0] << ", " << cand.local[1] << std::endl;
                // std::cout << "\tsurface: " << cand.surface_link << std::endl; 

                // Surface on which this measurement is
                detray::tracking_surface meas_surface(det, cand.surface_link);

                // std::cout << "\tsensitive? " << meas_surface.is_sensitive() << std::endl;
                // std::cout << "\thas material? " << meas_surface.has_material() << std::endl;
                
                // std::cout << "\thas material? " << meas_surface.has_material() << ", material: " << meas_surface.material_parameters(cand.local) << std::endl;
            }

            // Make a fitter state
            fitter.make_fitter(input_states);

            fitter.make_triplets();

            // Run fitter
            // fitter.fit(seed_param, fitter_state);

            // output_states.push_back(
                // std::move(fitter_state.m_fit_res),
                // std::move(fitter_state.m_fit_actor_state.m_track_states));

            // std::cout << "fitted chi2: " << output_states[i].header.chi2 << std::endl;

        }

        return output_states;
    }

    /// Config object
    config_type m_cfg;
};

}  // namespace traccc
