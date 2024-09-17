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
#include <string>

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

        // Matrix types
        using matrix_operator = detray::dmatrix_operator<algebra_type>;
        using size_type = detray::dsize_type<algebra_type>;
        template <size_type ROWS, size_type COLS>
        using matrix_type = detray::dmatrix<algebra_type, ROWS, COLS>;

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
            vecmem::vector<scalar> m_h_thet;
            vecmem::vector<scalar> m_h_phi;

            // Measurements for getting variances (hit shifts)
            // and surface orientation (scattering estimation)
            std::array<measurement, 3u> m_meas;
            
            // Estimated values
            scalar m_sigma_MS; // MS uncertainty
            scalar m_theta; // Polar angle 

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

                point2 loc_2d_0{meas_0.local[0], meas_0.local[1]};
                point2 loc_2d_1{meas_1.local[0], meas_1.local[1]};
                point2 loc_2d_2{meas_2.local[0], meas_2.local[1]};

                // Convert to global
                point3 glob_3d_0 = meas_0_sf.bound_to_global({}, loc_2d_0, {});
                point3 glob_3d_1 = meas_1_sf.bound_to_global({}, loc_2d_1, {});
                point3 glob_3d_2 = meas_2_sf.bound_to_global({}, loc_2d_2, {});

                std::cout << glob_3d_0[0] << " " << glob_3d_0[1] << " " << glob_3d_0[2] << std::endl;

                // Store triplet
                triplet t(glob_3d_0, glob_3d_1, glob_3d_2);
                
                // measurements copied here
                t.m_meas[0] = meas_0;
                t.m_meas[1] = meas_1;
                t.m_meas[2] = meas_2;
                
                // copy again
                m_triplets.push_back(t);

            }

            std::cout << m_triplets.size() << " triplets made\n";

        }


        /// Helper function - Linearize triplet
        ///
        /// Calculates triplet parameters by linearizing around circle solution
        /// @param t Triplet to linearize
        TRACCC_HOST_DEVICE void linearize_triplet(triplet& t) {

            // Curvature of circle in transverse plane

            std::cout << "Linearization:\n";

            vector3 x_01 {t.m_hit_1 - t.m_hit_0};
            vector3 x_12 {t.m_hit_2 - t.m_hit_1};

            scalar d_01 = getter::perp(x_01);
            scalar d_12 = getter::perp(x_12);
            scalar d_02 = getter::perp(vector3{t.m_hit_2 - t.m_hit_0});

            scalar z_01 = x_01[2];
            scalar z_12 = x_12[2];

            // TODO: x-prod evaluates to -ve, might have to be reversed
            scalar c_perp = 2.f * math::fabs((vector::cross(x_01, x_12))[2]) / (d_01 * d_12 * d_02); 

            std::cout << "\tc_perp " << c_perp << std::endl;


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

            std::cout << "\tc_3D_1C " << c_3D_1C << " c_3D_2C " << c_3D_2C << std::endl;

            // Polar angles
            scalar theta_1C = std::acos(z_01 * c_3D_1C / phi_1C);
            scalar theta_2C = std::acos(z_12 * c_3D_2C / phi_2C);
            t.m_theta = 0.5f * (theta_1C + theta_2C); // estimate polar angle

            vector2 tangent2D_norm = sin(t.m_theta) / getter::norm(tangent2D) * tangent2D;
            // track tangent normalized to 1
            vector3 tangent3D{tangent2D_norm[0], tangent2D_norm[1], static_cast<scalar>(cos(t.m_theta))};

            // Estimate MS-uncertainty
            
            detray::tracking_surface scat_sf(m_detector, t.m_meas[1].surface_link);

            // effective thickness
            scalar t_eff = mat_scatter / scat_sf.cos_angle({}, tangent3D, t.m_meas[1].local);

            auto scattering_unc = [](scalar curvature_3D, scalar eff_thickness, vector3 field_strength_vector) {
                return math::fabs(curvature_3D) * 45.f * math::sqrt(eff_thickness) * unit<scalar>::T / field_strength_vector[2] * (1.f + 0.038f * math::log(eff_thickness));
            };

            t.m_sigma_MS = scattering_unc((0.5f *(c_3D_1C + c_3D_2C)), t_eff, m_field.at(t.m_hit_1[0], t.m_hit_1[1], t.m_hit_1[2]));
            std::cout << "\tsigma_MS " << t.m_sigma_MS << std::endl;

            
            // Index parameters
            scalar n_1C = (d_01*d_01 + 4.f * z_01*z_01 * sin2_0p5_phi1C / phi2_1C) / (0.25f * d_01*d_01 * phi_1C * sin(phi_1C)/sin2_0p5_phi1C + 4.f * z_01*z_01 * sin2_0p5_phi1C/phi2_1C);
            scalar n_2C = (d_12*d_12 + 4.f * z_12*z_12 * sin2_0p5_phi2C / phi2_2C) / (0.25f * d_12*d_12 * phi_2C * sin(phi_2C)/sin2_0p5_phi2C + 4.f * z_12*z_12 * sin2_0p5_phi2C/phi2_2C);

            auto cot = [](scalar angle) { return cos(angle)/sin(angle); };
            // Triplet parameters
            t.m_phi_0 = 0.5f * (n_1C * phi_1C + n_2C * phi_2C);
            t.m_theta_0 = theta_2C - theta_1C + ((1.f -  n_2C) * cot(theta_2C) - (1.f - n_1C) * cot(theta_1C));
            t.m_rho_phi = -0.5f * (phi_1C * n_1C/c_3D_1C + phi_2C * n_2C/c_3D_2C);
            t.m_rho_theta = (1.f - n_1C) * cot(theta_1C) / c_3D_1C - (1.f - n_2C) * cot(theta_2C) / c_3D_2C;

            std::cout << "\tphi0 " << t.m_phi_0 << "  theta0 " << t.m_theta_0 << std::endl;
            std::cout << "\trho_phi " << t.m_rho_phi << "  rho_theta " << t.m_rho_theta << std::endl;
        }

        /// Helper function - Quick Linearize
        ///
        /// Faster calculation of theta_0 & phi_0
        /// using straight line trajectories 
        ///
        /// @param pos0 Position of hit0 (in global frame)
        /// @param pos1 Hit1
        /// @param pos2 Hit2
        TRACCC_HOST_DEVICE void quick_linearize(const vector3& pos0, const vector3& pos1, const vector3& pos2, scalar& phi_0, scalar& theta_0) {

            // make 2D vector from X, Y components of 3D vector
            auto perp_comp = [](const vector3& vec){ return vector2{vec[0], vec[1]}; };
            
            // Z-component magnitude of x-prod of 2D vectors
            auto cross_2d_z = [](const vector2& v1, const vector2& v2){ return math::fabs(v1[0]*v2[1] - v2[0]*v1[1]); };

            vector2 x_01 = perp_comp(pos1 - pos0);
            vector2 x_12 = perp_comp(pos2 - pos1);

            scalar d_01 = getter::norm(x_01);
            scalar d_12 = getter::norm(x_12);

            phi_0 = std::asin(cross_2d_z(x_01, x_12) / (d_01 * d_12));
            // std::cout << "Quick linearize:" << std::endl;
            // std::cout << "\tphi_0 " << phi_0 << std::endl;

            vector2 x_0_L{pos0[2], 0.f};
            vector2 x_1_L{pos1[2], d_01};
            vector2 x_2_L{pos2[2], d_01 + d_12};

            vector2 x_01_L = x_1_L - x_0_L;
            vector2 x_12_L = x_2_L - x_1_L;

            theta_0 = std::asin(cross_2d_z(x_01_L, x_12_L) / (getter::norm(x_01_L) * getter::norm(x_12_L)));
            // std::cout << "\ttheta_0 " << theta_0 << std::endl;
        }

        /// Helper function - Hit Position Derivatives
        ///
        /// Calulation of directional derivatives of
        /// triplet kinks w.r.t hit position shifts
        ///
        /// @param t Triplet
        TRACCC_HOST_DEVICE void calculate_pos_derivs(triplet& t) {

            // hits shifted by multiplier * sigma
            // in every hit uncertainty direction
            // scalar multiplier = 1.f;

            // need the position uncertainty directions for hits

            std::cout << "Hit position derivatives:\n";

            auto print_vec = [](vector3& v){ std::cout << v[0] << ", " << v[1] << ", " << v[2]; };

            auto R = [](scalar x, scalar y) { return math::sqrt(x*x + y*y); };

            detray::tracking_surface scat_sf(m_detector, t.m_meas[1].surface_link);
            context ctx;

            auto c_mass = scat_sf.centroid();

            auto c_geom = scat_sf.center(ctx);

            auto meas_loc = t.m_meas[1].local;
            auto meas_var = t.m_meas[1].variance;

            scalar phi_0_before = t.m_phi_0;
            scalar theta_0_before = t.m_theta_0;
            scalar phi_0_after;
            scalar theta_0_after;

            // Reserve space for derivative containers
            constexpr size_t max_dims = 3u;
            t.m_h_phi.reserve(3u * max_dims); // hits * max dims/hit
            t.m_h_thet.reserve(3u * max_dims);


            std::array<vector3, 3u> global_positions{t.m_hit_0, t.m_hit_1, t.m_hit_2};
            
            // Loop over measurements in triplet
            for (unsigned hit = 0; hit < 3; ++hit) {

                vector2 pos_loc = t.m_meas[hit].local;
                vector2 var_loc = t.m_meas[hit].variance;

                // Surface
                detray::tracking_surface sf(m_detector, t.m_meas[hit].surface_link);

                // over dimensions
                for (unsigned i = 0; i < max_dims; ++i) {

                    // Default derivative 0 for dimensions
                    // which don't exist for this measurement
                    if (i >= t.m_meas[hit].meas_dim) {
                        t.m_h_phi.push_back(0.f);
                        t.m_h_thet.push_back(0.f);
                        continue;
                    }

                    scalar sigma_i = math::sqrt(var_loc[i]);
                    
                    vector2 pos_shifted_loc = pos_loc;

                    // Shift (by the sigma in that direction)
                    pos_shifted_loc[i] = pos_loc[i] + sigma_i;

                    // In global frame
                    vector3 pos_shifted_glob = sf.bound_to_global({}, pos_shifted_loc, {}); 

                    global_positions[hit] = pos_shifted_glob;

                    // Get parameters with shifted hit
                    quick_linearize(global_positions[0], global_positions[1], global_positions[2], phi_0_after, theta_0_after);

                    t.m_h_phi.push_back((phi_0_after - phi_0_before) / sigma_i);
                    t.m_h_thet.push_back((theta_0_after - theta_0_before) / sigma_i);
                }

            }


            std::cout << "\tH_theta: ";
            for (unsigned j = 0; j < t.m_h_thet.size(); ++j) {
                std::cout << " " << t.m_h_thet[j];
            }
            std::cout << "\n\tH_phi: ";
            for (unsigned j = 0; j < t.m_h_phi.size(); ++j) {
                std::cout << " " << t.m_h_phi[j];
            }


            // vector3 vtx = scat_sf.local_vertices(0u);
            auto min_max_vtx = scat_sf.local_min_bounds();

            std::cout << "\nMiddle hit of triplet:\n";

            std::cout << "\tmeas: " << meas_loc[0] << ", " << meas_loc[1] << " +- " << meas_var[0] << ", " << meas_var[1] << "\n";

            std::cout << "\tmeas dimension: " << t.m_meas[1].meas_dim << "\n";

            std::cout << "\tsurface centroid: "; print_vec(c_mass); std::cout << "\n";

            std::cout << "\tsurface center: "; print_vec(c_geom); std::cout << " R: " << R(c_geom[0], c_geom[1]) << "\n";
            
            // std::cout << "vertex 0: "; print_vec(vtx); std ::cout << "\n";
            std::cout << "\tloc min bounds: " << min_max_vtx << std::endl;

            std::cout << "\tsurface details: " << scat_sf << std::endl;


        }

        /// Helper function - Global Fit
        ///
        /// Global fit of triplets
        ///
        TRACCC_HOST_DEVICE
        void do_global_fit() {

            // Allocate matrices with max possible sizes
            constexpr size_t max_dims = 3u;
            constexpr size_t max_nhits = 20u;
            constexpr size_t max_ntrips = max_nhits - 2u;
            constexpr size_t max_ndirs = max_dims * max_nhits; // max_dims = 3, better 2 ?

            // Actual number in this track
            const size_t N_triplets = m_triplets.size();
            const size_t N_hits = m_track_states.size();
            assert(N_hits <= max_nhits);
            assert(N_triplets == N_hits - 2u);


            // Make matrices/vectors

            // Triplet parameter vectors
            matrix_type<2u * max_ntrips, 1u> rho = matrix_operator().template zero<2u * max_ntrips, 1u>();
            matrix_type<2u * max_ntrips, 1u> psi = matrix_operator().template zero<2u * max_ntrips, 1u>();
            
            // Scattering & hit precision matrices
            matrix_type<2u * max_ntrips, 2u * max_ntrips> D_MS = matrix_operator().template identity<2u * max_ntrips, 2u * max_ntrips>();
            matrix_type<max_ndirs, max_ndirs> D_hit = matrix_operator().template zero<max_ndirs, max_ndirs>();

            // Position derivative matrices
            // matrix_type<max_ntrips, max_ndirs> H_theta = matrix_operator().template zero<max_ntrips, max_ndirs>();
            // matrix_type<max_ntrips, max_ndirs> H_phi = matrix_operator().template zero<max_ntrips, max_ndirs>();
            matrix_type<2u*max_ntrips, max_ndirs> H = matrix_operator().template zero<2u*max_ntrips, max_ndirs>();
            
            std::cout << getter::element(rho, 0u, 0u) << std::endl;

            
            // Fill matrices/vectors
            std::cout << "Filling matrices/vectors: \n";

            for (size_t i = 0; i < N_triplets; ++i) {
                std::cout << " Triplet " << i << std::endl;

                const triplet& t_i = m_triplets[i];

                getter::element(rho, i, 0u) = t_i.m_rho_theta;
                getter::element(rho, i + N_triplets, 0u) = t_i.m_rho_phi;

                getter::element(psi, i, 0u) = t_i.m_theta_0;
                getter::element(psi, i + N_triplets, 0u) = t_i.m_phi_0;

                scalar sigma2_MS = t_i.m_sigma_MS * t_i.m_sigma_MS;
                scalar sin2_theta = math::sin(t_i.m_theta);
                sin2_theta *= sin2_theta;
                getter::element(D_MS, i, i) = 1.f / sigma2_MS;
                getter::element(D_MS, i + N_triplets, i + N_triplets) = sin2_theta / sigma2_MS;


                // TODO: if max_dims is reduced to 2,
                // the manual assignements for the 3rd dimension
                // would go out of range !! -> loop over dimensions
                // might be unavoidable

                // (after unrolling loop over hits)
                // 1st Hit in triplet
                getter::element(H, i, i) = t_i.m_h_thet[0u];
                getter::element(H, i, max_nhits + i) = t_i.m_h_thet[1u];
                getter::element(H, i, 2u*max_nhits + i) = t_i.m_h_thet[2u];
                // 2nd Hit
                getter::element(H, i, i + 1u) = t_i.m_h_thet[max_dims*1u];
                getter::element(H, i, max_nhits + i + 1u) = t_i.m_h_thet[max_dims*1u + 1u];
                getter::element(H, i, 2u*max_nhits + i + 1u) = t_i.m_h_thet[max_dims*1u + 2u];
                // 3rd Hit
                getter::element(H, i, i + 2u) = t_i.m_h_thet[max_dims*2u];
                getter::element(H, i, max_nhits + i + 2u) = t_i.m_h_thet[max_dims*2u + 1u];
                getter::element(H, i, 2u*max_nhits + i + 2u) = t_i.m_h_thet[max_dims*2u + 2u];

                // 1st Hit
                getter::element(H, i + N_triplets, i) = t_i.m_h_phi[0u];
                getter::element(H, i + N_triplets, max_nhits + i) = t_i.m_h_phi[1u];
                getter::element(H, i + N_triplets, 2u*max_nhits + i) = t_i.m_h_phi[2u];
                // 2nd Hit
                getter::element(H, i + N_triplets, i + 1u) = t_i.m_h_phi[max_dims*1u];
                getter::element(H, i + N_triplets, max_nhits + i + 1u) = t_i.m_h_phi[max_dims*1u + 1u];
                getter::element(H, i + N_triplets, 2u*max_nhits + i + 1u) = t_i.m_h_phi[max_dims*1u + 2u];
                // 3rd Hit
                getter::element(H, i + N_triplets, i + 2u) = t_i.m_h_phi[max_dims*2u];
                getter::element(H, i + N_triplets, max_nhits + i + 2u) = t_i.m_h_phi[max_dims*2u + 1u];
                getter::element(H, i + N_triplets, 2u*max_nhits + i + 2u) = t_i.m_h_phi[max_dims*2u + 2u];
            

                // 1st Hit in triplet
                getter::element(D_hit, i * max_dims, i * max_dims) = 1.f / t_i.m_meas[0u].variance[0u];
                getter::element(D_hit, i * max_dims + 1u, i * max_dims + 1u) = 1.f / t_i.m_meas[0u].variance[1u];
                getter::element(D_hit, i * max_dims + 2u, i * max_dims + 2u) = 1.f; // dim does not exist

                // Only use the other two hits
                // for the last triplet to prevent
                // reassigning elements in marix
                if (i == N_triplets - 1u) {
                    // 2nd hit
                    getter::element(D_hit, (i + 1u) * max_dims, (i + 1u) * max_dims) = 1.f / t_i.m_meas[1u].variance[0u];
                    getter::element(D_hit, (i + 1u) * max_dims + 1u, (i + 1u) * max_dims + 1u) = 1.f / t_i.m_meas[1u].variance[1u];
                    getter::element(D_hit, (i + 1u) * max_dims + 2u, (i + 1u) * max_dims + 2u) = 1.f;

                    // 3rd hit
                    getter::element(D_hit, (i + 2u) * max_dims, (i + 2u) * max_dims) = 1.f / t_i.m_meas[2u].variance[0u];
                    getter::element(D_hit, (i + 2u) * max_dims + 1u, (i + 2u) * max_dims + 1u) = 1.f / t_i.m_meas[2u].variance[1u];
                    getter::element(D_hit, (i + 2u) * max_dims + 2u, (i + 2u) * max_dims + 2u) = 1.f;
                }

            } // done filling


            // Invert D_MS & D_hit diagonal matrices
            // (make diagonal elements corresponding to
            // unused 'objects' 1, to allow inversion)

            matrix_type<2u*max_ntrips, 2u*max_ntrips> D_MS_inv = matrix_operator().template identity<2u*max_ntrips, 2u*max_ntrips>();
            for (size_t j = 0u; j < 2u*max_ntrips; ++j) {

                // D_MS initialized as identity matrix
                getter::element(D_MS_inv, j, j) = 1.f / getter::element(D_MS, j, j); 
            }

            matrix_type<max_ndirs, max_ndirs> D_hit_inv = matrix_operator().template zero<max_ndirs, max_ndirs>();

            for (size_t j = 0u; j < max_ndirs; ++j) {
                
                // Leave 0 elements as it is
                // float comparison can be improved ?
                if (getter::element(D_hit, j , j) == 0.f) continue;
                getter::element(D_hit_inv, j, j) = 1.f / getter::element(D_hit, j, j);
            }

            // Triplet precision matrix
            // Note: diagonal elements in K_inv are 1
            // corresponding to unused 'objects', since
            // the same is true in D_MS_inv and those in
            // H * D_hit^-1 * H^T are 0

            matrix_type<2u*max_ntrips, 2u*max_ntrips> K_inv = D_MS_inv + H * D_hit_inv * matrix_operator().transpose(H);

            std::cout << "K_inv:\n";
            for (size_t r = 0u; r < 2u*max_ntrips; ++r) {
                for (size_t c = 0u; c < 2u*max_ntrips; ++c) {
                    std::cout << std::setw(12);
                    std::cout << getter::element(K_inv, r, c) << " ";
                }
                std::cout << std::endl;
            }

            // Matrix inversion
            matrix_type<2u*max_ntrips, 2u*max_ntrips> K = matrix_operator().inverse(K_inv);    

            matrix_type<1u, 1u> num = -1.f * matrix_operator().transpose(rho) * K * psi;
            matrix_type<1u, 1u> den = matrix_operator().transpose(rho) * K * rho;
            matrix_type<1u, 1u> psiT_K_psi = matrix_operator().transpose(psi) * K * psi;


            scalar c_3D = getter::element(num, 0u, 0u) / getter::element(den, 0u, 0u);

            scalar sigma_c_3D = 1.f / math::sqrt(getter::element(den, 0u, 0u));

            scalar chi2 = getter::element(psiT_K_psi, 0u, 0u) - (c_3D * c_3D) / (sigma_c_3D * sigma_c_3D); 

            std::cout << "Global fit: c_3D " << c_3D << "  sigma_c_3D " << sigma_c_3D << "  chi2 " << chi2 << std::endl;

            
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

            unsigned triplet_idx = 0;

            for (triplet& t : m_triplets) {
                std::cout << "Triplet " << triplet_idx << "\n";
                linearize_triplet(t);
                calculate_pos_derivs(t);
                ++triplet_idx;
                // break; // just one for debugging
            }

            do_global_fit();

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