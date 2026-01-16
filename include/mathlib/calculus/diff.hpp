#pragma once
#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "mathlib/core/error.hpp"

namespace mathlib::calculus {

	// Central difference derivative (good default)
	template <typename F, typename T>
	T derivative_central(F f, T x, T h = static_cast<T>(1e-6)) {
		static_assert(std::is_floating_point_v<T>, "derivative_central: T must be floating point");
		if (h <= T{}) throw core::domain_error("derivative_central(): h must be > 0");
		return (f(x + h) - f(x - h)) / (static_cast<T>(2) * h);
	}

	// Forward difference derivative (simpler, less accurate)
	template <typename F, typename T>
	T derivative_forward(F f, T x, T h = static_cast<T>(1e-6)) {
		static_assert(std::is_floating_point_v<T>, "derivative_forward: T must be floating point");
		if (h <= T{}) throw core::domain_error("derivative_forward(): h must be > 0");
		return (f(x + h) - f(x)) / h;
	}

	// Convenience default
	template <typename F, typename T>
	T derivative(F f, T x) {
		return derivative_central<F, T>(f, x);
	}

} // namespace mathlib::calculus
