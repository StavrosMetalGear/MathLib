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


    return 0;
}
