#pragma once
#include <cstddef>
#include <cmath>
#include <utility>

#include "mathlib/linalg/vector.hpp"
#include "mathlib/core/error.hpp"

namespace mathlib::calculus {

    template <std::size_t N, typename T>
    using Vec = mathlib::linalg::Vector<N, T>;

    // Simpson step for vector output
    template <typename F, std::size_t N, typename T>
    Vec<N, T> simpson_step_vec(F f, T a, T b) {
        const T c = (a + b) / static_cast<T>(2);
        return (b - a) / static_cast<T>(6) * (f(a) + static_cast<T>(4) * f(c) + f(b));
    }

    // Fixed-interval Simpson for vector output
    template <typename F, std::size_t N, typename T>
    Vec<N, T> integrate_simpson_vec(F f, T a, T b, std::size_t n = 1000) {
        if (n < 2) throw core::domain_error("integrate_simpson_vec(): n must be >= 2");
        if (n % 2 != 0) ++n;
        if (a == b) return Vec<N, T>{};
        if (b < a) std::swap(a, b);

        const T h = (b - a) / static_cast<T>(n);
        Vec<N, T> s = f(a) + f(b);

        for (std::size_t i = 1; i < n; ++i) {
            T x = a + static_cast<T>(i) * h;
            s = s + (i % 2 == 0 ? static_cast<T>(2) : static_cast<T>(4)) * f(x);
        }
        return s * (h / static_cast<T>(3));
    }

    // Adaptive Simpson for vector output
    template <typename F, std::size_t N, typename T>
    Vec<N, T> integrate_adaptive_simpson_vec(F f, T a, T b,
        T eps = static_cast<T>(1e-10),
        std::size_t max_recursion = 20) {
        if (eps <= T{}) throw core::domain_error("integrate_adaptive_simpson_vec(): eps must be > 0");
        if (a == b) return Vec<N, T>{};
        if (b < a) std::swap(a, b);

        const Vec<N, T> whole = simpson_step_vec<F, N, T>(f, a, b);

        struct Helper {
            static Vec<N, T> rec(F f, T a, T b, T eps, Vec<N, T> whole, std::size_t depth) {
                const T c = (a + b) / static_cast<T>(2);
                const Vec<N, T> left = simpson_step_vec<F, N, T>(f, a, c);
                const Vec<N, T> right = simpson_step_vec<F, N, T>(f, c, b);

                const Vec<N, T> delta = (left + right) - whole;

                // component-wise stopping: use max abs component
                T max_abs = T{};
                for (std::size_t i = 0; i < N; ++i) max_abs = std::max(max_abs, std::abs(delta[i]));

                if (depth == 0 || max_abs <= static_cast<T>(15) * eps) {
                    return left + right + delta / static_cast<T>(15);
                }

                return rec(f, a, c, eps / static_cast<T>(2), left, depth - 1) +
                    rec(f, c, b, eps / static_cast<T>(2), right, depth - 1);
            }
        };

        return Helper::rec(f, a, b, eps, whole, max_recursion);
    }

} // namespace mathlib::calculus
#pragma once
