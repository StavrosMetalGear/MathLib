#pragma once
#include <cstddef>
#include <cmath>
#include <limits>

#include "mathlib/core/error.hpp"
#include "mathlib/calculus/diff.hpp"

namespace mathlib::calculus {

    // Bisection: requires f(a) and f(b) have opposite signs.
    template <typename F, typename T>
    T root_bisection(F f, T a, T b,
        T eps = static_cast<T>(1e-12),
        std::size_t max_iter = 200) {
        if (eps <= T{}) throw core::domain_error("root_bisection(): eps must be > 0");
        T fa = f(a), fb = f(b);
        if (fa == T{}) return a;
        if (fb == T{}) return b;
        if ((fa > T{} && fb > T{}) || (fa < T{} && fb < T{})) {
            throw core::domain_error("root_bisection(): f(a) and f(b) must have opposite signs");
        }

        for (std::size_t it = 0; it < max_iter; ++it) {
            T m = (a + b) / static_cast<T>(2);
            T fm = f(m);

            if (std::abs(fm) <= eps || (b - a) / static_cast<T>(2) <= eps) return m;

            if ((fa > T{} && fm > T{}) || (fa < T{} && fm < T{})) {
                a = m; fa = fm;
            }
            else {
                b = m; fb = fm;
            }
        }
        return (a + b) / static_cast<T>(2);
    }

    // Newton: fast but needs decent initial guess; uses numeric derivative by default.
    template <typename F, typename T>
    T root_newton(F f, T x0,
        T eps = static_cast<T>(1e-12),
        std::size_t max_iter = 50,
        T h = static_cast<T>(1e-6)) {
        if (eps <= T{}) throw core::domain_error("root_newton(): eps must be > 0");

        T x = x0;
        for (std::size_t it = 0; it < max_iter; ++it) {
            T fx = f(x);
            if (std::abs(fx) <= eps) return x;

            T dfx = derivative_central<F, T>(f, x, h);
            if (dfx == T{}) throw core::domain_error("root_newton(): derivative is zero");

            T step = fx / dfx;
            x = x - step;

            if (std::abs(step) <= eps) return x;
        }
        return x;
    }

} // namespace mathlib::calculus
