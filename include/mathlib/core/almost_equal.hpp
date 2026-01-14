#pragma once
#include <cmath>
#include <type_traits>
#include <algorithm>

namespace mathlib::core {

    // Robust float compare: abs + relative tolerance
    template <typename T>
    constexpr bool almost_equal(T a, T b,
        T rel_tol = static_cast<T>(1e-12),
        T abs_tol = static_cast<T>(1e-12)) {
        static_assert(std::is_floating_point_v<T>, "almost_equal requires floating point");

        // Handle infinities exactly
        if (a == b) return true;

        const T diff = std::abs(a - b);
        const T norm = std::max(std::abs(a), std::abs(b));
        return diff <= std::max(abs_tol, rel_tol * norm);
    }

} // namespace mathlib::core
#pragma once
