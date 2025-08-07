#include "vector2d.h"
#include <cmath>

Vector2D::Vector2D(double x, double y) : x(x), y(y) {}

Vector2D Vector2D::operator+(const Vector2D& other) const {
    return Vector2D(x + other.x, y + other.y);
}

Vector2D Vector2D::operator-(const Vector2D& other) const {
    return Vector2D(x - other.x, y - other.y);
}

double Vector2D::dot(const Vector2D& other) const {
    return x * other.x + y * other.y;
}

double Vector2D::magnitude() const {
    return std::sqrt(x * x + y * y);
}
