#include <gtest/gtest.h>
#include <cmath>

#include "mathlib/core/almost_equal.hpp"
#include "mathlib/core/constants.hpp"
#include "mathlib/linalg/vector.hpp"
#include "mathlib/calculus/integrate_vec.hpp"
#include "mathlib/calculus/grad.hpp"
#include "mathlib/calculus/root.hpp"
#include "mathlib/ode/solvers.hpp"

TEST(VectorValuedIntegration, IntegrateCosSin) {
    using mathlib::core::almost_equal;
    using mathlib::linalg::Vector;

    auto f = [](double x) {
        return Vector<2, double>{ std::cos(x), std::sin(x) };
        };

    const double pi = mathlib::core::pi_v<double>;
    auto I = mathlib::calculus::integrate_adaptive_simpson_vec<decltype(f), 2, double>(f, 0.0, pi, 1e-10);

    // ∫0^pi cos(x) dx = 0
    // ∫0^pi sin(x) dx = 2
    EXPECT_TRUE(almost_equal(I[0], 0.0, 1e-9, 1e-9));
    EXPECT_TRUE(almost_equal(I[1], 2.0, 1e-9, 1e-9));
}

TEST(Gradient, Quadratic) {
    using mathlib::core::almost_equal;
    using mathlib::linalg::Vector;

    // f(x,y) = x^2 + 3y^2 -> grad = (2x, 6y)
    auto f = [](const Vector<2, double>& v) {
        return v[0] * v[0] + 3.0 * v[1] * v[1];
        };

    Vector<2, double> x{ 2.0, -1.0 };
    auto g = mathlib::calculus::gradient<decltype(f), 2, double>(f, x);

    EXPECT_TRUE(almost_equal(g[0], 4.0, 1e-6, 1e-6));
    EXPECT_TRUE(almost_equal(g[1], -6.0, 1e-6, 1e-6));
}

TEST(RootFinding, BisectionAndNewton) {
    using mathlib::core::almost_equal;

    // root of x^2 - 2 = 0 is sqrt(2)
    auto f = [](double x) { return x * x - 2.0; };

    double r1 = mathlib::calculus::root_bisection<decltype(f), double>(f, 1.0, 2.0, 1e-12);
    double r2 = mathlib::calculus::root_newton<decltype(f), double>(f, 1.5, 1e-12);

    EXPECT_TRUE(almost_equal(r1, std::sqrt(2.0), 1e-10, 1e-10));
    EXPECT_TRUE(almost_equal(r2, std::sqrt(2.0), 1e-10, 1e-10));
}

TEST(ODE, RK4ExpGrowth) {
    using mathlib::core::almost_equal;

    // y' = y, y(0)=1 => y(t)=e^t
    auto f = [](double /*t*/, double y) { return y; };

    auto traj = mathlib::ode::solve_rk4<decltype(f), double, double>(f, 0.0, 1.0, 1.0, 1e-3);
    double y1 = traj.back().second;

    EXPECT_TRUE(almost_equal(y1, std::exp(1.0), 1e-3, 1e-3));
}

TEST(ODE, RK45ExpGrowth) {
    using mathlib::core::almost_equal;

    auto f = [](double /*t*/, double y) { return y; };
    auto traj = mathlib::ode::solve_rk45<decltype(f), double, double>(f, 0.0, 1.0, 1.0, 1e-2, 1e-9);
    double y1 = traj.back().second;

    EXPECT_TRUE(almost_equal(y1, std::exp(1.0), 1e-6, 1e-6));
}
