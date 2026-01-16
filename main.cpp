#include <iostream>
#include "vector2d.h"
#include <iostream>
#include "vector3d.h"

int main() {
    Vector3D a(1, 0, 0);
    Vector3D b(0, 1, 0);

    Vector3D c = a.cross(b);
    std::cout << "a x b = (" << c.x << ", " << c.y << ", " << c.z << ")" << std::endl;
    Vector3D v(3, 4, 0);
    Vector3D unit = v.normalize();
    std::cout << "Normalized vector = (" << unit.x << ", " << unit.y << ", " << unit.z << ")" << std::endl;
    mathlib::linalg::Matrix<2, 2, double> A2{ 2,1,5,7 };
    mathlib::linalg::Vector<2, double> b2{ 11,13 };
    auto x2 = mathlib::linalg::solve(A2, b2);
    std::cout << "solve(A,b) = (" << x2[0] << ", " << x2[1] << ")\n";
    std::cout << "d/dx sin(x) at 0 = "
        << mathlib::calculus::derivative([](double x) { return std::sin(x); }, 0.0)
        << "\n";

    std::cout << "Integral sin(x) from 0 to pi = "
        << mathlib::calculus::integrate_adaptive_simpson([](double x) { return std::sin(x); }, 0.0, mathlib::core::pi_v<double>)
        << "\n";


    return 0;
}
