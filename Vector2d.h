#pragma once

class Vector2D {
public:
    double x, y;

    Vector2D(double x = 0, double y = 0);

    Vector2D operator+(const Vector2D& other) const;
    Vector2D operator-(const Vector2D& other) const;
    double dot(const Vector2D& other) const;
    double magnitude() const;
};

