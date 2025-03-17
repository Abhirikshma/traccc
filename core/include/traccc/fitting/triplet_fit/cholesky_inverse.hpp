// Custom matrix inversion algorithm
// using Cholseky decomposition


// #include <iostream>

namespace traccc {


    template <typename algebra_t>
    class cholesky_inverse {

        public:
        // Detector type
        // using detector_type = typename navigator_t::detector_type;

        // Algebra type
        // using algebra_type = typename detector_type::algebra_type;
        using algebra_type = algebra_t;

        // Matrix types
        using matrix_operator = detray::dmatrix_operator<algebra_type>;
        using size_type = detray::dsize_type<algebra_type>;
        template <size_type ROWS, size_type COLS>
        using matrix_type = detray::dmatrix<algebra_type, ROWS, COLS>;


        template <size_type N>
        TRACCC_HOST_DEVICE inline matrix_type<N, N> operator()(
            const matrix_type<N, N>& m) const {

            // std::cout << "N " << N << std::endl;

            // Inverse matrix
            matrix_type<N, N> inv = matrix_operator().template zero<N, N>();

            // L (lower traingular)
            matrix_type<N, N> L = matrix_operator().template zero<N, N>();

            // std::cout << "finding L\n";
            // Cholesky–Banachiewicz algorithm
            for (size_type i = 0; i < N; ++i) {
                for (size_type j = 0; j <= i; ++j) {

                    scalar sum = 0.f;

                    for (size_type k = 0; k < j; ++k) {
                        sum += getter::element(L, i, k) * getter::element(L, j, k);
                    }

                    if (i == j)
                    getter::element(L, i, j) = math::sqrt(getter::element(m, i, i) - sum);

                    else
                    getter::element(L, i, j) = 1.f / getter::element(L, j, j) * (getter::element(m, i, j) - sum);
                }
            }

            // L^T (upper triangular)
            matrix_type<N, N> L_T = matrix_operator().transpose(L);

            // Y (LY = I) - solved by forward substitution
            matrix_type<N, N> Y = matrix_operator().template zero<N, N>();

            // std::cout << "forward subs\n";

            for (size_type j = 0; j < N; ++j) {
                for (size_type i = 0; i < N; ++i) {
                    
                    scalar sum = 0.f;

                    for (size_type k = 0; k < i; ++k)
                        sum += getter::element(L, i, k) * getter::element(Y, k, j);


                    getter::element(Y, i, j) = 1.f / getter::element(L, i, i) * ((i == j ? 1.f : 0.f) - sum);
                }
            }

            // X (L^T X = Y) - solved by backward substitution
            // X is the inverse, i.e. inv

            // std::cout << "backward subs\n";

            for (size_type j = N; j-- > 0; ) {
                for (size_type i = N; i-- > 0; ) {

                    // std::cout << "j " << j << " i " << i << std::endl;
                    
                    scalar sum = 0.f;

                    for (size_type k = N; k-- > i; ) {
                        // std::cout << "k " << k << std::endl;
                        sum += getter::element(L_T, i, k) * getter::element(inv, k, j);
                    }

                    getter::element(inv, i, j) = 1.f / getter::element(L_T, i, i) * (getter::element(Y, i, j) - sum);
                }
            }
            

            // std::cout << "inversion using cholesky\n";

            return inv;
        }


    };
}