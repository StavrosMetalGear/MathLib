#pragma once
#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>

namespace mathlib::linalg {

    template <std::size_t R, std::size_t C, typename T = double>
    struct Matrix {
        static_assert(R > 0 && C > 0, "Matrix dimensions must be > 0");
        static_assert(std::is_arithmetic_v<T>, "Matrix<T>: T must be arithmetic");

        // row-major storage
        std::array<T, R* C> a{};

        constexpr Matrix() = default;

        // Allow: Matrix<2,3>{ 1,2,3,4,5,6 }
        constexpr Matrix(std::initializer_list<T> init) {
            if (init.size() != R * C) {
                throw std::invalid_argument("Matrix initializer_list size mismatch");
            }
            std::size_t i = 0;
            for (auto& x : init) a[i++] = x;
        }

        constexpr T& operator()(std::size_t r, std::size_t c) { return a[r * C + c]; }
        constexpr const T& operator()(std::size_t r, std::size_t c) const { return a[r * C + c]; }

        static constexpr Matrix identity() requires (R == C) {
            Matrix I;
            for (std::size_t i = 0; i < R; ++i) I(i, i) = T{ 1 };
            return I;
        }
    };

    // Matrix + Matrix
    template <std::size_t R, std::size_t C, typename T>
    constexpr Matrix<R, C, T> operator+(const Matrix<R, C, T>& A, const Matrix<R, C, T>& B) {
        Matrix<R, C, T> out;
        for (std::size_t i = 0; i < R * C; ++i) out.a[i] = A.a[i] + B.a[i];
        return out;
    }

    template <std::size_t R, std::size_t C, typename T>
    constexpr Matrix<R, C, T> operator-(const Matrix<R, C, T>& A, const Matrix<R, C, T>& B) {
        Matrix<R, C, T> out;
        for (std::size_t i = 0; i < R * C; ++i) out.a[i] = A.a[i] - B.a[i];
        return out;
    }

    // Matrix * Matrix
    template <std::size_t R, std::size_t K, std::size_t C, typename T>
    constexpr Matrix<R, C, T> operator*(const Matrix<R, K, T>& A, const Matrix<K, C, T>& B) {
        Matrix<R, C, T> out;
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t c = 0; c < C; ++c) {
                T sum{};
                for (std::size_t k = 0; k < K; ++k) sum += A(r, k) * B(k, c);
                out(r, c) = sum;
            }
        }
        return out;
    }

    // Transpose
    template <std::size_t R, std::size_t C, typename T>
    constexpr Matrix<C, R, T> transpose(const Matrix<R, C, T>& A) {
        Matrix<C, R, T> out;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t c = 0; c < C; ++c)
                out(c, r) = A(r, c);
        return out;
    }

} // namespace mathlib::linalg
#pragma once
