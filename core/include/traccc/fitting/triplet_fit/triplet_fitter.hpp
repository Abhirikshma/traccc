/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_candidate.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state.hpp"
#include "traccc/fitting/fitting_config.hpp"

// detray include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"

// System include(s).
#include <limits>
#include <iostream>

namespace traccc {

    /// Triplet fitting algorithm to fit a single track

    // we dont need the stepper/navigator now but
    // the detector and bfield types are obtained from these,
    // as in the original Kalman filter code
    template <typename stepper_t, typename navigator_t>
    class triplet_fitter {

        public:
        // Detector type
        using detector_type = typename navigator_t::detector_type;

        // Algebra type
        using algebra_type = typename detector_type::algebra_type;

        // Vector type
        template <typename T>
        using vector_type = typename navigator_t::template vector_type<T>;

        // Configuration type
        using config_type = fitting_config;

        // Field type
        using bfield_type = typename stepper_t::magnetic_field_type;

        /// Constructor with a detector
        ///
        /// @param det the detector object
        ///
        TRACCC_HOST_DEVICE
        triplet_fitter(const detector_type& det, const bfield_type& field, const config_type& cfg)
                    : m_detector(det), m_field(field), m_cfg(cfg) {}

        // Triplet struct
        struct triplet {

            // positions of three hits
            point3 hit_0;
            point3 hit_1;
            point3 hit_2;

            // triplet parameters
            scalar phi_0;
            scalar theta_0;
            scalar rho_phi;
            scalar rho_theta;

            // hit position derivatives
            vector_type<algebra_type> h_thet;
            vector_type<algebra_type> h_phi;

            // extra intermediate information if needed

        };

        /// Helper function - Initialize fitter
        TRACCC_HOST_DEVICE
        void make_fitter(const vector_type<track_state<algebra_type>>& in_track_states) {

            m_track_states = in_track_states;
        }

        /// Helper function - Make triplets 
        ///
        /// Makes triplets from consecutive measurements on the track
        ///
        TRACCC_HOST_DEVICE 
        void make_triplets() {
            
            // loop over measurements (track states) in candidate
            for (const track_state<algebra_type>& state : m_track_states) {
                
                measurement meas = state.get_measurement();
                std::cout << "measurement local pos: " << meas.local[0] << ", " << meas.local[1] << std::endl;

                // point3 local_3d{meas.local[0], meas.local[1], 0.f};

                // Get surface
                // detray::tracking_surface meas_sf(m_detector, meas.surface_link);

                // Convert to global
                

            }

            // get rid of unused variable error!
            triplet t;
            m_triplets.push_back(t);
            

        }

        /// State -> better name ?
        /// Can remove this in favour of having the variables contained here
        /// as members of the enclosing triplet_fitter class ?
        /// 
        /// Contains:
        /// (1) container of triplets of the track being fitted
        /// (2) final fitting result
        /// 
        // struct state {

        //     /// Constructor makes triplets
        //     ///
        //     /// @param in_track_states vector of input track states
        //     ///
        //     TRACCC_HOST_DEVICE
        //     state(const vector_type<track_state<algebra_type>> && in_track_states, const triplet_fitter& fitter) 
        //     : m_track_states(std::move(in_track_states)) {
        //         make_triplets(std::move(in_track_states), fitter.m_detector, this->m_triplets);
        //     }

        //     /// Over-loaded constructor
        //     ///
        //     /// FIXME!
        //     TRACCC_HOST_DEVICE
        //     state(const vector_type<track_state<algebra_type>> && in_track_states)
        //     : m_track_states(std::move(in_track_states)) {}

        //     /// Store track states (just in case)
        //     vector_type<track_state<algebra_type>> m_track_states;
            
        //     /// Vector of triplets
        //     vector_type<triplet> m_triplets;

        //     /// Final fitting result per track
        //     fitting_result<algebra_type> m_fit_res;
        // };

        /// Helper function - Linearize triplet
        ///
        /// Calculates triplet parameters by linearizing around circle solution
        ///
        TRACCC_HOST_DEVICE void linearize_triplet(triplet& t) {

        }
        
        /// Run the fitter
        /// 
        /// 
        TRACCC_HOST_DEVICE void fit() {

            // main function, calls all other functions

            // make triplets -> can be a separate function
            // e.g. fitter_state for Kalman

            // calculate triplet parameters

            // calculate position derivatives

            // global fit

        }
        
        /// Store track states
        vector_type<track_state<algebra_type>> m_track_states;
        
        private:
        // Detector object
        const detector_type& m_detector;
        // Field object
        const bfield_type m_field;
        
        /// Vector of triplets
        vector_type<triplet> m_triplets;
        

        // Configuration obeject
        config_type m_cfg;
    
    };

} // namespace traccc