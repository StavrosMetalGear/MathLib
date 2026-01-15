#pragma once
#include <cstddef>
#include <stdexcept>
#include <cmath>

#include "mathlib/linalg/matrix.hpp"
#include "mathlib/linalg/vector.hpp"
#include "mathlib/core/almost_equal.hpp"
#include "mathlib/core/error.hpp"

namespace mathlib::linalg {

    // Solve A x = b using Gaussian elimination with partial pivoting.
    // Works for fixed-size square matrices Matrix<N,N,T> and Vector<N,T>.
    template <std::size_t N, typename T>
    Vector<N, T> solve(Matrix<N, N, T> A, Vector<N, T> b,
        T pivot_eps = static_cast<T>(1e-12)) {
        static_assert(N > 0, "solve<N>: N must be > 0");

        // Forward elimination
        for (std::size_t k = 0; k < N; ++k) {
            // Find pivot row p with max |A(p,k)| for p>=k
            std::size_t pivot = k;
            T max_abs = std::abs(A(k, k));
            for (std::size_t i = k + 1; i < N; ++i) {
                T v = std::abs(A(i, k));
                if (v > max_abs) {
                    max_abs = v;
                    pivot = i;
                }
            }

            // Check near-singular pivot
            if (max_abs <= pivot_eps) {
                throw core::domain_error("solve(): matrix is singular or ill-conditioned (pivot ~ 0)");
            }

            // Swap pivot row into place
            if (pivot != k) {
                for (std::size_t j = k; j < N; ++j) {
                    std::swap(A(k, j), A(pivot, j));
                }
                std::swap(b[k], b[pivot]);
            }

            // Eliminate below pivot
            for (std::size_t i = k + 1; i < N; ++i) {
                T factor = A(i, k) / A(k, k);
                A(i, k) = T{}; // exactly zero out
                for (std::size_t j = k + 1; j < N; ++j) {
                    A(i, j) -= factor * A(k, j);
                }
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        Vector<N, T> x;
        for (std::size_t i = N; i-- > 0;) {
            T sum = b[i];
            for (std::size_t j = i + 1; j < N; ++j) {
                sum -= A(i, j) * x[j];
            }
            if (std::abs(A(i, i)) <= pivot_eps) {
                throw core::domain_error("solve(): matrix is singular or ill-conditioned (diag ~ 0)");
            }
            x[i] = sum / A(i, i);
        }
        return x;
    }

    // Helper: compute A*x (useful for tests and examples)
    template <std::size_t N, typename T>
    Vector<N, T> mul(const Matrix<N, N, T>& A, const Vector<N, T>& x) {
        Vector<N, T> y;
        for (std::size_t r = 0; r < N; ++r) {
            T sum{};
            for (std::size_t c = 0; c < N; ++c) sum += A(r, c) * x[c];
            y[r] = sum;
        }
        return y;
    }

} // namespace mathlib::linalg
#pragma once
