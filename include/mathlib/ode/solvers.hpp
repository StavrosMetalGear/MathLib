#pragma once
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

#include "mathlib/core/error.hpp"

namespace mathlib::ode {

    // A trajectory is a list of (t, y) samples.
    template <typename State, typename T>
    using Trajectory = std::vector<std::pair<T, State>>;

    // Utility: max-norm for error control (works for double and your Vector)
    inline double max_norm(double x) { return std::abs(x); }

    template <std::size_t N, typename T>
    T max_norm(const mathlib::linalg::Vector<N, T>& v) {
        T m{};
        for (std::size_t i = 0; i < N; ++i) m = std::max(m, std::abs(v[i]));
        return m;
    }

    // Euler (fixed step)
    template <typename F, typename State, typename T>
    Trajectory<State, T> solve_euler(F f, T t0, State y0, T t1, T h) {
        if (h <= T{}) throw mathlib::core::domain_error("solve_euler(): h must be > 0");
        if (t1 < t0) std::swap(t0, t1);

        Trajectory<State, T> out;
        out.reserve(static_cast<std::size_t>((t1 - t0) / h) + 2);

        T t = t0;
        State y = y0;
        out.push_back({ t, y });

        while (t < t1) {
            T step = std::min(h, t1 - t);
            y = y + step * f(t, y);
            t += step;
            out.push_back({ t, y });
        }
        return out;
    }

    // RK4 (fixed step)
    template <typename F, typename State, typename T>
    Trajectory<State, T> solve_rk4(F f, T t0, State y0, T t1, T h) {
        if (h <= T{}) throw mathlib::core::domain_error("solve_rk4(): h must be > 0");
        if (t1 < t0) std::swap(t0, t1);

        Trajectory<State, T> out;
        out.reserve(static_cast<std::size_t>((t1 - t0) / h) + 2);

        T t = t0;
        State y = y0;
        out.push_back({ t, y });

        while (t < t1) {
            T step = std::min(h, t1 - t);

            auto k1 = f(t, y);
            auto k2 = f(t + step / 2, y + (step / 2) * k1);
            auto k3 = f(t + step / 2, y + (step / 2) * k2);
            auto k4 = f(t + step, y + step * k3);

            y = y + (step / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
            t += step;
            out.push_back({ t, y });
        }

        return out;
    }

    // RK45 adaptive (Dormand–Prince 5(4))
    template <typename F, typename State, typename T>
    Trajectory<State, T> solve_rk45(F f, T t0, State y0, T t1,
        T h0 = static_cast<T>(1e-2),
        T eps = static_cast<T>(1e-9),
        T h_min = static_cast<T>(1e-10),
        T h_max = static_cast<T>(1.0)) {
        if (h0 <= T{}) throw mathlib::core::domain_error("solve_rk45(): h0 must be > 0");
        if (eps <= T{}) throw mathlib::core::domain_error("solve_rk45(): eps must be > 0");
        if (t1 < t0) std::swap(t0, t1);

        Trajectory<State, T> out;
        out.reserve(1024);

        T t = t0;
        State y = y0;
        T h = std::clamp(h0, h_min, h_max);

        out.push_back({ t, y });

        while (t < t1) {
            if (h < h_min) throw mathlib::core::domain_error("solve_rk45(): step underflow (h < h_min)");
            if (t + h > t1) h = t1 - t;

            // Dormand–Prince coefficients
            auto k1 = f(t, y);
            auto k2 = f(t + h * (T(1) / 5), y + h * (T(1) / 5) * k1);
            auto k3 = f(t + h * (T(3) / 10), y + h * (T(3) / 40) * k1 + h * (T(9) / 40) * k2);
            auto k4 = f(t + h * (T(4) / 5), y + h * (T(44) / 45) * k1 + h * (T(-56) / 15) * k2 + h * (T(32) / 9) * k3);
            auto k5 = f(t + h * (T(8) / 9), y + h * (T(19372) / 6561) * k1 + h * (T(-25360) / 2187) * k2 + h * (T(64448) / 6561) * k3 + h * (T(-212) / 729) * k4);
            auto k6 = f(t + h, y + h * (T(9017) / 3168) * k1 + h * (T(-355) / 33) * k2 + h * (T(46732) / 5247) * k3 + h * (T(49) / 176) * k4 + h * (T(-5103) / 18656) * k5);

            // 5th order solution
            State y5 = y + h * (T(35) / 384) * k1 + h * (T(500) / 1113) * k3 + h * (T(125) / 192) * k4 + h * (T(-2187) / 6784) * k5 + h * (T(11) / 84) * k6;

            // 4th order solution (error estimate)
            auto k7 = f(t + h, y5);
            State y4 = y + h * (T(5179) / 57600) * k1 + h * (T(7571) / 16695) * k3 + h * (T(393) / 640) * k4
                + h * (T(-92097) / 339200) * k5 + h * (T(187) / 2100) * k6 + h * (T(1) / 40) * k7;

            State err = y5 - y4;
            T err_norm = static_cast<T>(max_norm(err));

            // Accept/reject step
            if (err_norm <= eps || err_norm == T{}) {
                t += h;
                y = y5;
                out.push_back({ t, y });
            }

            // Update h (classic controller)
            T safety = static_cast<T>(0.9);
            T power = static_cast<T>(0.2); // 1/5
            T factor = (err_norm == T{}) ? static_cast<T>(5) : safety * std::pow(eps / err_norm, power);
            factor = std::clamp(factor, static_cast<T>(0.2), static_cast<T>(5.0));
            h = std::clamp(h * factor, h_min, h_max);
        }

        return out;
    }

} // namespace mathlib::ode
