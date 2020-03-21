#pragma once

#include "wake/core/types.hpp"
#include "wake/core/constants.hpp"
#include <array>

namespace wake {

// D2Q9 Lattice Model
// 9 discrete velocities in 2D space:
// c_0 = (0, 0)
// c_1..4 = (1,0), (0,1), (-1,0), (0,-1)
// c_5..8 = (1,1), (-1,1), (-1,-1), (1,-1)
struct D2Q9 {
    static constexpr int Q = 9;
    static constexpr int D = 2;

    // Discrete velocity vectors [Q][D]
    static constexpr std::array<int_coord_t, Q> CX = {0,  1,  0, -1,  0,  1, -1, -1,  1};
    static constexpr std::array<int_coord_t, Q> CY = {0,  0,  1,  0, -1,  1,  1, -1, -1};

    // Lattice quadrature weights
    static constexpr std::array<real_t, Q> W = {
        4.0f / 9.0f,
        1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f
    };

    // Opposite direction indices for bounce-back boundary conditions
    static constexpr std::array<int, Q> OPPOSITE = {0, 3, 4, 1, 2, 7, 8, 5, 6};
};

// D3Q19 Lattice Model
// 19 discrete velocities in 3D space:
// Rest (1), Faces (6), Edges (12)
struct D3Q19 {
    static constexpr int Q = 19;
    static constexpr int D = 3;

    // Discrete velocity vectors
    static constexpr std::array<int_coord_t, Q> CX = {
        0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0
    };
    static constexpr std::array<int_coord_t, Q> CY = {
        0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1
    };
    static constexpr std::array<int_coord_t, Q> CZ = {
        0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1
    };

    // Lattice quadrature weights
    // w_0 = 1/3, w_1..6 = 1/18, w_7..18 = 1/36
    static constexpr std::array<real_t, Q> W = {
        1.0f / 3.0f,
        1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f
    };

    // Opposite direction indices for bounce-back
    static constexpr std::array<int, Q> OPPOSITE = {
        0,  2,  1,  4,  3,  6,  5,  10, 9,  8,  7,  14, 13, 12, 11, 18, 17, 16, 15
    };
};

// D3Q27 Lattice Model
// 27 discrete velocities in 3D space:
// Rest (1), Faces (6), Edges (12), Corners (8)
struct D3Q27 {
    static constexpr int Q = 27;
    static constexpr int D = 3;

    static constexpr std::array<int_coord_t, Q> CX = {
        0,
         1, -1,  0,  0,  0,  0,
         1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,
         1, -1,  1, -1,  1, -1,  1, -1
    };
    static constexpr std::array<int_coord_t, Q> CY = {
        0,
         0,  0,  1, -1,  0,  0,
         1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,
         1,  1, -1, -1,  1,  1, -1, -1
    };
    static constexpr std::array<int_coord_t, Q> CZ = {
        0,
         0,  0,  0,  0,  1, -1,
         0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
         1,  1,  1,  1, -1, -1, -1, -1
    };

    // Weights: w_0 = 8/27, w_1..6 = 2/27, w_7..18 = 1/54, w_19..26 = 1/216
    static constexpr std::array<real_t, Q> W = {
        8.0f / 27.0f,
        2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f, 2.0f / 27.0f,
        1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f,
        1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f,
        1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f, 1.0f / 54.0f,
        1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f,
        1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f, 1.0f / 216.0f
    };

    static constexpr std::array<int, Q> OPPOSITE = {
        0,
        2,  1,  4,  3,  6,  5,
        10, 9,  8,  7,  14, 13, 12, 11, 18, 17, 16, 15,
        22, 21, 20, 19, 26, 25, 24, 23
    };
};

} // namespace wake
