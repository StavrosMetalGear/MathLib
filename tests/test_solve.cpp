#include <gtest/gtest.h>

#include "mathlib/linalg/solve.hpp"
#include "mathlib/core/almost_equal.hpp"

TEST(LinalgSolve, Solves2x2) {
    using namespace mathlib::linalg;

    Matrix<2, 2, double> A{
      2, 1,
      5, 7
    };
    Vector<2, double> b{ 11, 13 };

    // Expected solution: x = 64/9, y = -29/9
    auto x = solve(A, b);

    EXPECT_TRUE(mathlib::core::almost_equal(x[0], 64.0 / 9.0));
    EXPECT_TRUE(mathlib::core::almost_equal(x[1], -29.0 / 9.0));
}

TEST(LinalgSolve, Solves3x3AndChecksResidual) {
    using namespace mathlib::linalg;
    using mathlib::core::almost_equal;

    Matrix<3, 3, double> A{
      3, 2, -1,
      2, -2, 4,
      -1, 0.5, -1
    };
    Vector<3, double> b{ 1, -2, 0 };

    auto x = solve(A, b);
    auto r = mul(A, x);

    EXPECT_TRUE(almost_equal(r[0], b[0], 1e-12, 1e-12));
    EXPECT_TRUE(almost_equal(r[1], b[1], 1e-12, 1e-12));
    EXPECT_TRUE(almost_equal(r[2], b[2], 1e-12, 1e-12));
}

TEST(LinalgSolve, SingularThrows) {
    using namespace mathlib::linalg;

    // Rows are linearly dependent -> singular
    Matrix<2, 2, double> A{
      1, 2,
      2, 4
    };
    Vector<2, double> b{ 3, 6 };

    EXPECT_THROW((void)solve(A, b), mathlib::core::domain_error);
}
