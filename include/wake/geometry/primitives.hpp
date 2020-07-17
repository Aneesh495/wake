#pragma once

#include "wake/core/types.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/grid/lattice_grid_3d.hpp"

namespace wake {

class GeometryPrimitives {
public:
    // 2D Primitives
    static void rasterize_circle(LatticeGrid2D& grid, real_t cx, real_t cy, real_t radius);
    static void rasterize_rectangle(LatticeGrid2D& grid, real_t x0, real_t y0, real_t x1, real_t y1);
    static void rasterize_rotated_plate(LatticeGrid2D& grid, real_t cx, real_t cy, real_t length, real_t thickness, real_t angle_deg);
    static void rasterize_porous_array(LatticeGrid2D& grid, real_t start_x, real_t end_x, real_t spacing, real_t obstacle_radius);

    // 3D Primitives
    static void voxelize_sphere(LatticeGrid3D& grid, real_t cx, real_t cy, real_t cz, real_t radius);
    static void voxelize_cylinder_z(LatticeGrid3D& grid, real_t cx, real_t cy, real_t radius);
    static void voxelize_box(LatticeGrid3D& grid, real_t x0, real_t y0, real_t z0, real_t x1, real_t y1, real_t z1);
    static void voxelize_ellipsoid(LatticeGrid3D& grid, real_t cx, real_t cy, real_t cz, real_t rx, real_t ry, real_t rz);
};

} // namespace wake
