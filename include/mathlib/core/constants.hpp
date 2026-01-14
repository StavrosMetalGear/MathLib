#pragma once
#include <limits>
#include <type_traits>

namespace mathlib::core {

	template <typename T>
	inline constexpr T pi_v = static_cast<T>(3.141592653589793238462643383279502884L);

	template <typename T>
	inline constexpr T e_v = static_cast<T>(2.718281828459045235360287471352662498L);

	template <typename T>
	constexpr T epsilon() {
		static_assert(std::is_floating_point_v<T>, "epsilon<T>() requires floating point");
		return std::numeric_limits<T>::epsilon();
	}

} // namespace mathlib::core
#pragma once
