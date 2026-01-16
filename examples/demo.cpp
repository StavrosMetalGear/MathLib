#include <iostream>
#include "mathlib/linalg/vector.hpp"
#include "mathlib/linalg/matrix.hpp"
#include "mathlib/linalg/solve.hpp"
#include "mathlib/calculus/diff.hpp"
#include "mathlib/calculus/integrate.hpp"


int main() {
    using namespace mathlib::linalg;

    Vec3 a{ 1, 0, 0 };
    Vec3 b{ 0, 1, 0 };

    auto c = cross(a, b);
    std::cout << "cross(a,b) = (" << c[0] << ", " << c[1] << ", " << c[2] << ")\n";

    Matrix<2, 2> A{ 1, 2,
                   3, 4 };
    Matrix<2, 2> B{ 2, 0,
                   1, 2 };
    auto C = A * B;

    std::cout << "A*B =\n";
    std::cout << C(0, 0) << " " << C(0, 1) << "\n";
    std::cout << C(1, 0) << " " << C(1, 1) << "\n";
}
