#include <gtest/gtest.h>
#include "mathlib/linalg/vector.hpp"
#include "mathlib/linalg/matrix.hpp"

TEST(Vector, DotAndNorm) {
	using mathlib::linalg::Vector;
	Vector<3> v{ 3, 4, 0 };
	EXPECT_DOUBLE_EQ(v.norm2(), 25.0);
	EXPECT_DOUBLE_EQ(v.norm(), 5.0);
}

TEST(Vector, Cross) {
	using namespace mathlib::linalg;
	Vec3 a{ 1,0,0 };
	Vec3 b{ 0,1,0 };
	auto c = cross(a, b);
	EXPECT_DOUBLE_EQ(c[0], 0.0);
	EXPECT_DOUBLE_EQ(c[1], 0.0);
	EXPECT_DOUBLE_EQ(c[2], 1.0);
}

TEST(Matrix, Multiply2x2) {
	using mathlib::linalg::Matrix;
	Matrix<2, 2> A{ 1,2,3,4 };
	Matrix<2, 2> B{ 2,0,1,2 };
	auto C = A * B;

	EXPECT_DOUBLE_EQ(C(0, 0), 1 * 2 + 2 * 1);
	EXPECT_DOUBLE_EQ(C(0, 1), 1 * 0 + 2 * 2);
	EXPECT_DOUBLE_EQ(C(1, 0), 3 * 2 + 4 * 1);
	EXPECT_DOUBLE_EQ(C(1, 1), 3 * 0 + 4 * 2);
}
