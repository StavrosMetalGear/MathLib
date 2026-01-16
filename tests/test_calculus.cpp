#include <gtest/gtest.h>
#include <cmath>

#include "mathlib/calculus/diff.hpp"
#include "mathlib/calculus/integrate.hpp"
#include "mathlib/core/almost_equal.hpp"
#include "mathlib/core/constants.hpp"

TEST(Calculus, DerivativeSinAt0) {
	using mathlib::core::almost_equal;
	auto f = [](double x) { return std::sin(x); };

	double d = mathlib::calculus::derivative(f, 0.0);
	EXPECT_TRUE(almost_equal(d, 1.0, 1e-8, 1e-8));
}

TEST(Calculus, DerivativePolynomial) {
	using mathlib::core::almost_equal;
	auto f = [](double x) { return x * x * x; }; // x^3

	double d = mathlib::calculus::derivative(f, 2.0);
	EXPECT_TRUE(almost_equal(d, 12.0, 1e-6, 1e-6)); // 3x^2 = 12
}

TEST(Calculus, IntegrateSimpsonSin0Pi) {
	using mathlib::core::almost_equal;
	const double pi = mathlib::core::pi_v<double>;
	auto f = [](double x) { return std::sin(x); };

	double I = mathlib::calculus::integrate_simpson<decltype(f), double>(f, 0.0, pi, 2000);
	EXPECT_TRUE(almost_equal(I, 2.0, 1e-8, 1e-8));
}

TEST(Calculus, IntegrateAdaptiveSimpsonExp0_1) {
	using mathlib::core::almost_equal;
	auto f = [](double x) { return std::exp(x); };

	// ∫_0^1 e^x dx = e - 1
	double I = mathlib::calculus::integrate_adaptive_simpson<decltype(f), double>(f, 0.0, 1.0, 1e-12);
	EXPECT_TRUE(almost_equal(I, std::exp(1.0) - 1.0, 1e-10, 1e-10));
}
