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
        /// @param field magnetic field
        /// @param cfg fitter configuration
        TRACCC_HOST_DEVICE
        triplet_fitter(const detector_type& det, const bfield_type& field, const config_type& cfg)
                    : m_detector(det), m_field(field), m_cfg(cfg) {}

        // Triplet struct
        struct triplet {

            /// Default construct
            triplet() = default;

            /// Construct with three hits
            ///
            /// @param hit_0 position of first hit
            /// @param hit_1 position of second hit
            /// @param hit_2 position of third hit
            triplet(const point3& hit_0, const point3& hit_1, const point3& hit_2)
                : m_hit_0(hit_0), m_hit_1(hit_1), m_hit_2(hit_2) { }

            // positions of three hits
            point3 m_hit_0;
            point3 m_hit_1;
            point3 m_hit_2;

            // triplet parameters
            scalar m_phi_0;
            scalar m_theta_0;
            scalar m_rho_phi;
            scalar m_rho_theta;

            // hit position derivatives
            vector_type<algebra_type> m_h_thet;
            vector_type<algebra_type> m_h_phi;

            // extra intermediate information if needed
            measurement m_meas_1; // measurement for middle hit
            scalar m_sigma_MS; // estimation of MS uncertainty

        };

        /// Helper function - Initialize fitter
        ///
        /// @param in_track_states input track states (from measurements)
        TRACCC_HOST_DEVICE
        void init_fitter(const vector_type<track_state<algebra_type>>& in_track_states) {

            m_track_states = in_track_states;

            m_triplets.clear();
        }

        /// Helper function - Make triplets 
        ///
        /// Makes triplets from consecutive measurements on the track
        ///
        TRACCC_HOST_DEVICE 
        void make_triplets() {
            
            std::size_t n_triplets = m_track_states.size() - 2;

            m_triplets.reserve(n_triplets);

            // loop over measurements (track states) in candidate
            for (std::size_t i = 0; i < n_triplets; ++i) {

                std::cout << "\t" << i << "-th triplet\n";
                
                // Get track states (and measurements)
                const track_state<algebra_type>& state_0 = m_track_states[i];
                const track_state<algebra_type>& state_1 = m_track_states[i+1];
                const track_state<algebra_type>& state_2 = m_track_states[i+2];
                
                const measurement& meas_0 = state_0.get_measurement();
                const measurement& meas_1 = state_1.get_measurement();
                const measurement& meas_2 = state_2.get_measurement();

                // Get surfaces
                detray::tracking_surface meas_0_sf(m_detector, meas_0.surface_link);
                detray::tracking_surface meas_1_sf(m_detector, meas_1.surface_link);
                detray::tracking_surface meas_2_sf(m_detector, meas_2.surface_link);

                // Convert to global
                context ctx;
                vector3 dir{};

                point2 loc_2d_0{meas_0.local[0], meas_0.local[1]};
                point2 loc_2d_1{meas_1.local[0], meas_1.local[1]};
                point2 loc_2d_2{meas_2.local[0], meas_2.local[1]};

                point3 glob_3d_0 = meas_0_sf.bound_to_global(ctx, loc_2d_0, dir);
                point3 glob_3d_1 = meas_1_sf.bound_to_global(ctx, loc_2d_1, dir);
                point3 glob_3d_2 = meas_2_sf.bound_to_global(ctx, loc_2d_2, dir);

                // Store triplet
                triplet t(glob_3d_0, glob_3d_1, glob_3d_2);
                t.m_meas_1 = meas_1;
                m_triplets.push_back(t);

            }

            std::cout << m_triplets.size() << " triplets made\n";

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
        /// @param t Triplet to linearize
        TRACCC_HOST_DEVICE void linearize_triplet(triplet& t) {

            // Curvature of circle in transverse plane

            vector3 x_01 {t.m_hit_1 - t.m_hit_0};
            vector3 x_12 {t.m_hit_2 - t.m_hit_1};

            scalar d_01 = getter::perp(x_01);
            scalar d_12 = getter::perp(x_12);
            scalar d_02 = getter::perp(vector3{t.m_hit_2 - t.m_hit_0});

            scalar z_01 = x_01[2];
            scalar z_12 = x_12[2];

            // TODO: x-prod evaluates to -ve, might have to be reversed
            scalar c_perp = 2.f * math::fabs((vector::cross(x_01, x_12))[2]) / (d_01 * d_12 * d_02); 

            std::cout << "c_perp " << c_perp << std::endl;


            // Direction of track at scattering plane (using hits 0 & 2)
            vector2 m{0.5f * (t.m_hit_0[0] + t.m_hit_2[0]), 0.5f * (t.m_hit_0[1] + t.m_hit_2[1])};
            
            vector2 n{(t.m_hit_2[1] - t.m_hit_0[1]) / d_02, (t.m_hit_0[0] - t.m_hit_2[0]) / d_02};

            scalar perp_d = math::sqrt(1.f / (c_perp * c_perp) - 0.25f * (d_02 * d_02));

            // two centres possible
            std::array<vector2, 2u> c;

            c[0] = vector2{m[0] + n[0] * perp_d, m[1] + n[1] * perp_d};
            c[1] = vector2{m[0] - n[0] * perp_d, m[1] - n[1] * perp_d};
            
            vector2 x1{t.m_hit_1[0], t.m_hit_1[1]};
            
            assert(getter::norm(x1 - m) != 0.f); // 3 hits must not lie on a straight line

            // choose the correct one
            vector2 c_correct{0.f, 0.f};
            for (const vector2& c_i : c) {
                if (vector::dot(x1 - m, c_i - m) < 0.f) {
                    c_correct = c_i;
                    break;
                }
            }
            assert(getter::norm(c_correct) > 0.f); // no correct centre of circle found
            vector2 r1 = x1 - c_correct;
            vector2 tangent2D{r1[1], -1.f * r1[0]};

            // tangent direction along trajectory
            if (vector::dot(tangent2D, vector2{x_12[0], x_12[1]}) < 0.f)
                tangent2D = -1.f * tangent2D;



            // Parameters of the arc segments
            
            // Azimuthal (bending) angles
            scalar phi_1C = 2.f * std::asin(0.5f * d_01 * c_perp);
            scalar phi_2C = 2.f * std::asin(0.5f * d_12 * c_perp);
            
            scalar phi2_1C = phi_1C * phi_1C;
            scalar phi2_2C = phi_2C * phi_2C;
            
            scalar sin2_0p5_phi1C = sin(0.5f * phi_1C);
            sin2_0p5_phi1C *= sin2_0p5_phi1C;
            scalar sin2_0p5_phi2C = sin(0.5f * phi_2C);
            sin2_0p5_phi2C *= sin2_0p5_phi2C;

            // 3D curvatures
            scalar c_3D_1C = phi_1C / math::sqrt(z_01*z_01 + 0.25f * d_01*d_01 * phi2_1C / sin2_0p5_phi1C);
            scalar c_3D_2C = phi_2C / math::sqrt(z_12*z_12 + 0.25f * d_12*d_12 * phi2_2C / sin2_0p5_phi2C);

            std::cout << "c_3D_1C " << c_3D_1C << " c_3D_2C " << c_3D_2C << std::endl;

            // Polar angles
            scalar theta_1C = std::acos(z_01 * c_3D_1C / phi_1C);
            scalar theta_2C = std::acos(z_12 * c_3D_2C / phi_2C);
            scalar theta_est = 0.5f * (theta_1C + theta_2C); // estimate polar angle

            vector2 tangent2D_norm = sin(theta_est) / getter::norm(tangent2D) * tangent2D;
            // track tangent normalized to 1
            vector3 tangent3D{tangent2D_norm[0], tangent2D_norm[1], static_cast<scalar>(cos(theta_est))};

            // Estimate MS-uncertainty
            context ctx;
            detray::tracking_surface scat_sf(m_detector, t.m_meas_1.surface_link);

            // effective thickness
            scalar t_eff = mat_scatter / scat_sf.cos_angle(ctx, tangent3D, t.m_meas_1.local);

            auto scattering_unc = [](scalar curvature_3D, scalar eff_thickness, vector3 field_strength_vector) {
                return math::fabs(curvature_3D) * 45.f * math::sqrt(eff_thickness) * unit<scalar>::T / field_strength_vector[2] * (1.f + 0.038f * math::log(eff_thickness));
            };

            t.m_sigma_MS = scattering_unc((0.5f *(c_3D_1C + c_3D_2C)), t_eff, m_field.at(t.m_hit_1[0], t.m_hit_1[1], t.m_hit_1[2]));
            std::cout << "sigma_MS " << t.m_sigma_MS << std::endl;

            
            // Index parameters
            scalar n_1C = (d_01*d_01 + 4.f * z_01*z_01 * sin2_0p5_phi1C / phi2_1C) / (0.25f * d_01*d_01 * phi_1C * sin(phi_1C)/sin2_0p5_phi1C + 4.f * z_01*z_01 * sin2_0p5_phi1C/phi2_1C);
            scalar n_2C = (d_12*d_12 + 4.f * z_12*z_12 * sin2_0p5_phi2C / phi2_2C) / (0.25f * d_12*d_12 * phi_2C * sin(phi_2C)/sin2_0p5_phi2C + 4.f * z_12*z_12 * sin2_0p5_phi2C/phi2_2C);

            auto cot = [](scalar angle) { return cos(angle)/sin(angle); };
            // Triplet parameters
            t.m_phi_0 = 0.5f * (n_1C * phi_1C + n_2C * phi_2C);
            t.m_theta_0 = theta_2C - theta_1C + ((1.f -  n_2C) * cot(theta_2C) - (1.f - n_1C) * cot(theta_1C));
            t.m_rho_phi = -0.5f * (phi_1C * n_1C/c_3D_1C + phi_2C * n_2C/c_3D_2C);
            t.m_rho_theta = (1.f - n_1C) * cot(theta_1C) / c_3D_1C - (1.f - n_2C) * cot(theta_2C) / c_3D_2C;

            std::cout << "phi0 " << t.m_phi_0 << "  theta0 " << t.m_theta_0 << std::endl;
            std::cout << "rho_phi " << t.m_rho_phi << "  rho_theta " << t.m_rho_theta << std::endl;
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

            for (triplet& t : m_triplets) {
                linearize_triplet(t);
                break; // just one for debugging
            }

        }
        
        
        private:

        // Hard-coded material
        // TODO: get from surface after re-mapping
        scalar mat_scatter = 0.02f;

        // Detector context type
        using context = typename detector_type::geometry_context;

        // Detector object
        const detector_type& m_detector;
        // Field object
        const bfield_type m_field;
        
        /// Vector of triplets
        vector_type<triplet> m_triplets;
        // Track states
        vector_type<track_state<algebra_type>> m_track_states;
        
        // Configuration object
        config_type m_cfg;
    
    };

} // namespace traccc