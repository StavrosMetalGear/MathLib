#include <gtest/gtest.h>
#include "mathlib/core/constants.hpp"
#include "mathlib/core/almost_equal.hpp"

TEST(Core, Constants) {
	EXPECT_GT(mathlib::core::pi_v<double>, 3.14);
	EXPECT_LT(mathlib::core::pi_v<double>, 3.15);
}

TEST(Core, AlmostEqual) {
	using mathlib::core::almost_equal;

	EXPECT_TRUE(almost_equal(1.0, 1.0 + 1e-13));
	EXPECT_FALSE(almost_equal(1.0, 1.0 + 1e-6));
	EXPECT_TRUE(almost_equal(0.0, 1e-13, 1e-12, 1e-12));
}
