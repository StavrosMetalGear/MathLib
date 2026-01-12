#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>

namespace mathlib::linalg {

    template <std::size_t N, typename T = double>
    struct Vector {
        static_assert(N > 0, "Vector dimension N must be > 0");
        static_assert(std::is_arithmetic_v<T>, "Vector<T>: T must be arithmetic");

        std::array<T, N> v{};

        constexpr Vector() = default;

        // Allow: Vector<3>{1,2,3}
        constexpr Vector(std::initializer_list<T> init) {
            if (init.size() != N) {
                throw std::invalid_argument("Vector initializer_list size mismatch");
            }
            std::size_t i = 0;
            for (auto& x : init) v[i++] = x;
        }

        constexpr T& operator[](std::size_t i) { return v[i]; }
        constexpr const T& operator[](std::size_t i) const { return v[i]; }

        // Vector + Vector
        friend constexpr Vector operator+(const Vector& a, const Vector& b) {
            Vector out;
            for (std::size_t i = 0; i < N; ++i) out[i] = a[i] + b[i];
            return out;
        }

        friend constexpr Vector operator-(const Vector& a, const Vector& b) {
            Vector out;
            for (std::size_t i = 0; i < N; ++i) out[i] = a[i] - b[i];
            return out;
        }

        // Scalar ops
        friend constexpr Vector operator*(const Vector& a, T s) {
            Vector out;
            for (std::size_t i = 0; i < N; ++i) out[i] = a[i] * s;
            return out;
        }
        friend constexpr Vector operator*(T s, const Vector& a) { return a * s; }

        friend constexpr Vector operator/(const Vector& a, T s) {
            if (s == T{}) throw std::invalid_argument("Vector division by zero scalar");
            Vector out;
            for (std::size_t i = 0; i < N; ++i) out[i] = a[i] / s;
            return out;
        }

        // Dot product
        friend constexpr T dot(const Vector& a, const Vector& b) {
            T sum{};
            for (std::size_t i = 0; i < N; ++i) sum += a[i] * b[i];
            return sum;
        }

        // Squared length, length
        constexpr T norm2() const { return dot(*this, *this); }

        T norm() const {
            using std::sqrt;
            return sqrt(norm2());
        }

        Vector normalized(T eps = static_cast<T>(1e-12)) const {
            T n = norm();
            if (n <= eps) throw std::domain_error("Cannot normalize near-zero vector");
            return (*this) / n;
        }
    };

    // Common aliases
    using Vec2 = Vector<2, double>;
    using Vec3 = Vector<3, double>;
    using Vec4 = Vector<4, double>;

    // Cross product only for 3D
    template <typename T>
    constexpr Vector<3, T> cross(const Vector<3, T>& a, const Vector<3, T>& b) {
        return Vector<3, T>{
            a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
        };
    }

} // namespace mathlib::linalg
#pragma once
