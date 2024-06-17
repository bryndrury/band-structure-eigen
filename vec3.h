#pragma once

#include <vector>
#include <iostream>
#include <cmath>

struct vec3
{
    double x, y, z;

    vec3() : x(0), y(0), z(0) {}
    vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    vec3(vec3 const &v) : x(v.x), y(v.y), z(v.z) {}

    // Operators
    vec3 operator+(const vec3 &v) const
    {
        return vec3(x + v.x, y + v.y, z + v.z);
    }
    vec3 operator-(const vec3 &v) const
    {
        return vec3(x - v.x, y - v.y, z - v.z);
    }
    vec3 operator*(const double i) const
    {
        return vec3(x * i, y * i, z * i);
    }
    vec3 operator/(const double i) const
    {
        return vec3(x / i, y / i, z / i);
    }

    vec3 &operator=(const vec3 &v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }
    vec3 &operator+=(const vec3 &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    vec3 &operator-=(const vec3 &v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    vec3 &operator*=(const double i)
    {
        x *= i;
        y *= i;
        z *= i;
        return *this;
    }
    vec3 &operator/=(const double i)
    {
        x /= i;
        y /= i;
        z /= i;
        return *this;
    }

    vec3 cross(const vec3 &v) const
    {
        return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    double dot(const vec3 &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    double norm() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    void out() const
    {
        std::cout << x << " " << y << " " << z << std::endl;
    }
};