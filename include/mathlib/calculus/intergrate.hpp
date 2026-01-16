#pragma once
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "mathlib/core/error.hpp"

namespace mathlib::calculus {

    // Simpson's rule with even n subintervals
    template <typename F, typename T>
    T integrate_simpson(F f, T a, T b, std::size_t n = 1000) {
        static_assert(std::is_floating_point_v<T>, "integrate_simpson: T must be floating point");
        if (n < 2) throw core::domain_error("integrate_simpson(): n must be >= 2");
        if (n % 2 != 0) ++n; // make even
        if (a == b) return T{};
        if (b < a) std::swap(a, b);

        const T h = (b - a) / static_cast<T>(n);
        T s = f(a) + f(b);

        for (std::size_t i = 1; i < n; ++i) {
            T x = a + static_cast<T>(i) * h;
            s += (i % 2 == 0 ? static_cast<T>(2) : static_cast<T>(4)) * f(x);
        }
        return s * (h / static_cast<T>(3));
    }

    // Internal: one Simpson step on [a,b]
    template <typename F, typename T>
    T simpson_step(F f, T a, T b) {
        const T c = (a + b) / static_cast<T>(2);
        return (b - a) / static_cast<T>(6) * (f(a) + static_cast<T>(4) * f(c) + f(b));
    }

    // Adaptive Simpson's rule
    template <typename F, typename T>
    T integrate_adaptive_simpson(F f, T a, T b,
        T eps = static_cast<T>(1e-10),
        std::size_t max_recursion = 20) {
        static_assert(std::is_floating_point_v<T>, "integrate_adaptive_simpson: T must be floating point");
        if (eps <= T{}) throw core::domain_error("integrate_adaptive_simpson(): eps must be > 0");
        if (a == b) return T{};
        if (b < a) std::swap(a, b);

        const T whole = simpson_step<F, T>(f, a, b);

        // recursive lambda (C++17-friendly)
        struct Helper {
            static T rec(F f, T a, T b, T eps, T whole, std::size_t depth) {
                const T c = (a + b) / static_cast<T>(2);
                const T left = simpson_step<F, T>(f, a, c);
                const T right = simpson_step<F, T>(f, c, b);
                const T delta = left + right - whole;

                // if good enough or depth exhausted
                if (depth == 0 || std::abs(delta) <= static_cast<T>(15) * eps) {
                    // Richardson extrapolation correction
                    return left + right + delta / static_cast<T>(15);
                }
                return rec(f, a, c, eps / static_cast<T>(2), left, depth - 1) +
                    rec(f, c, b, eps / static_cast<T>(2), right, depth - 1);
            }
        };

        return Helper::rec(f, a, b, eps, whole, max_recursion);
    }

} // namespace mathlib::calculus
