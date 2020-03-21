#pragma once

#include "wake/core/types.hpp"

namespace wake {

// Fundamental Lattice Boltzmann dimensionless physical constants
// Lattice speed of sound: cs = 1 / sqrt(3)
inline constexpr real_t CS = 0.5773502691896258f;
inline constexpr real_t CS2 = 1.0f / 3.0f;
inline constexpr real_t INV_CS2 = 3.0f;
inline constexpr real_t INV_2_CS4 = 4.5f;
inline constexpr real_t INV_2_CS2 = 1.5f;

// Reference density for weakly compressible isothermal flow
inline constexpr real_t RHO_0 = 1.0f;

// SIMD vectorization alignment requirements (64 bytes for AVX-512 and cache line alignment)
inline constexpr std::size_t SIMD_ALIGNMENT = 64;

// Default Smagorinsky subgrid-scale turbulence constant
inline constexpr real_t DEFAULT_SMAGORINSKY_CONSTANT = 0.14f;

// Default WALE turbulence model constant
inline constexpr real_t DEFAULT_WALE_CONSTANT = 0.50f;

// Universal gravitational acceleration default (lattice units)
inline constexpr real_t DEFAULT_GRAVITY = 9.80665f;

// Math constants
inline constexpr real_t PI = 3.14159265358979323846f;
inline constexpr real_t TWO_PI = 6.28318530717958647692f;

} // namespace wake
