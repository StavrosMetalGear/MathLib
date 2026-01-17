#pragma once
#include <cstddef>
#include <type_traits>

#include "mathlib/linalg/vector.hpp"
#include "mathlib/core/error.hpp"

namespace mathlib::calculus {

    // Gradient of scalar function f: R^N -> R using central differences
    template <typename F, std::size_t N, typename T>
    mathlib::linalg::Vector<N, T> gradient(F f,
        const mathlib::linalg::Vector<N, T>& x,
        T h = static_cast<T>(1e-6)) {
        static_assert(std::is_floating_point_v<T>, "gradient: T must be floating point");
        if (h <= T{}) throw mathlib::core::domain_error("gradient(): h must be > 0");

        mathlib::linalg::Vector<N, T> g{};

        for (std::size_t i = 0; i < N; ++i) {
            auto xp = x;
            auto xm = x;
            xp[i] += h;
            xm[i] -= h;
            g[i] = (f(xp) - f(xm)) / (static_cast<T>(2) * h);
        }
        return g;
    }

} // namespace mathlib::calculus

