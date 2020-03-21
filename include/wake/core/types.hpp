#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <array>

namespace wake {

// Scalar floating point precision for host and solver coordinates
using real_t = float;

// Integer index and coordinate representations
using index_t = std::size_t;
using int_coord_t = int32_t;

// Vector representations in 2D and 3D
struct Vec2 {
    real_t x{0.0f};
    real_t y{0.0f};

    constexpr Vec2() = default;
    constexpr Vec2(real_t x_, real_t y_) : x(x_), y(y_) {}

    constexpr Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    constexpr Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    constexpr Vec2 operator*(real_t s) const { return {x * s, y * s}; }
    constexpr Vec2 operator/(real_t s) const { return {x / s, y / s}; }

    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
    Vec2& operator-=(const Vec2& o) { x -= o.x; y -= o.y; return *this; }
    Vec2& operator*=(real_t s) { x *= s; y *= s; return *this; }

    constexpr real_t dot(const Vec2& o) const { return x * o.x + y * o.y; }
    constexpr real_t norm_sq() const { return x * x + y * y; }
};

struct Vec3 {
    real_t x{0.0f};
    real_t y{0.0f};
    real_t z{0.0f};

    constexpr Vec3() = default;
    constexpr Vec3(real_t x_, real_t y_, real_t z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator*(real_t s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(real_t s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(real_t s) { x *= s; y *= s; z *= s; return *this; }

    constexpr real_t dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    constexpr real_t norm_sq() const { return x * x + y * y + z * z; }

    constexpr Vec3 cross(const Vec3& o) const {
        return {
            y * o.z - z * o.y,
            z * o.x - x * o.z,
            x * o.y - y * o.x
        };
    }
};

// Node boundary type identifiers
enum class CellType : uint8_t {
    Fluid = 0,
    Solid = 1,
    Inflow = 2,
    Outflow = 3,
    Periodic = 4,
    Symmetry = 5,
    MovingWall = 6
};

// Collision model options
enum class CollisionType : uint8_t {
    BGK = 0,
    TRT = 1,
    MRT = 2,
    Entropic = 3
};

// Turbulence modeling options
enum class TurbulenceModel : uint8_t {
    None = 0,
    Smagorinsky = 1,
    WALE = 2
};

// Memory layout options
enum class MemoryLayout : uint8_t {
    AoS = 0,     // Array of Structures
    SoA = 1,     // Structure of Arrays (SIMD preferred)
    AoSoA = 2    // Array of Structures of Arrays (Cache blocked)
};

} // namespace wake
