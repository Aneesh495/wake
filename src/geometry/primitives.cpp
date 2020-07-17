#include "wake/geometry/primitives.hpp"
#include <cmath>
#include <algorithm>

namespace wake {

void GeometryPrimitives::rasterize_circle(LatticeGrid2D& grid, real_t cx, real_t cy, real_t radius) {
    grid.set_solid_circle(cx, cy, radius);
}

void GeometryPrimitives::rasterize_rectangle(LatticeGrid2D& grid, real_t x0, real_t y0, real_t x1, real_t y1) {
    grid.set_solid_box(static_cast<int>(x0), static_cast<int>(y0), static_cast<int>(x1), static_cast<int>(y1));
}

void GeometryPrimitives::rasterize_rotated_plate(
    LatticeGrid2D& grid, real_t cx, real_t cy, real_t length, real_t thickness, real_t angle_deg
) {
    real_t rad = angle_deg * (PI / 180.0f);
    real_t cos_a = std::cos(rad);
    real_t sin_a = std::sin(rad);

    real_t half_len = 0.5f * length;
    real_t half_thick = 0.5f * thickness;

    int min_x = std::max(0, static_cast<int>(cx - half_len - thickness));
    int max_x = std::min(grid.nx() - 1, static_cast<int>(cx + half_len + thickness));
    int min_y = std::max(0, static_cast<int>(cy - half_len - thickness));
    int max_y = std::min(grid.ny() - 1, static_cast<int>(cy + half_len + thickness));

    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            real_t dx = static_cast<real_t>(x) - cx;
            real_t dy = static_cast<real_t>(y) - cy;

            // Rotate coordinates into plate local frame
            real_t local_x = dx * cos_a + dy * sin_a;
            real_t local_y = -dx * sin_a + dy * cos_a;

            if (std::abs(local_x) <= half_len && std::abs(local_y) <= half_thick) {
                grid.set_cell_type(x, y, CellType::Solid);
            }
        }
    }
}

void GeometryPrimitives::rasterize_porous_array(
    LatticeGrid2D& grid, real_t start_x, real_t end_x, real_t spacing, real_t obstacle_radius
) {
    for (real_t x = start_x; x <= end_x; x += spacing) {
        for (real_t y = spacing; y < static_cast<real_t>(grid.ny()) - spacing; y += spacing) {
            rasterize_circle(grid, x, y, obstacle_radius);
        }
    }
}

void GeometryPrimitives::voxelize_sphere(LatticeGrid3D& grid, real_t cx, real_t cy, real_t cz, real_t radius) {
    grid.set_solid_sphere(cx, cy, cz, radius);
}

void GeometryPrimitives::voxelize_cylinder_z(LatticeGrid3D& grid, real_t cx, real_t cy, real_t radius) {
    grid.set_solid_cylinder_z(cx, cy, radius);
}

void GeometryPrimitives::voxelize_box(
    LatticeGrid3D& grid, real_t x0, real_t y0, real_t z0, real_t x1, real_t y1, real_t z1
) {
    grid.set_solid_box(
        static_cast<int>(x0), static_cast<int>(y0), static_cast<int>(z0),
        static_cast<int>(x1), static_cast<int>(y1), static_cast<int>(z1)
    );
}

void GeometryPrimitives::voxelize_ellipsoid(
    LatticeGrid3D& grid, real_t cx, real_t cy, real_t cz, real_t rx, real_t ry, real_t rz
) {
    int min_x = std::max(0, static_cast<int>(cx - rx - 1));
    int max_x = std::min(grid.nx() - 1, static_cast<int>(cx + rx + 1));
    int min_y = std::max(0, static_cast<int>(cy - ry - 1));
    int max_y = std::min(grid.ny() - 1, static_cast<int>(cy + ry + 1));
    int min_z = std::max(0, static_cast<int>(cz - rz - 1));
    int max_z = std::min(grid.nz() - 1, static_cast<int>(cz + rz + 1));

    real_t inv_rx2 = 1.0f / (rx * rx);
    real_t inv_ry2 = 1.0f / (ry * ry);
    real_t inv_rz2 = 1.0f / (rz * rz);

    for (int z = min_z; z <= max_z; ++z) {
        for (int y = min_y; y <= max_y; ++y) {
            for (int x = min_x; x <= max_x; ++x) {
                real_t dx = static_cast<real_t>(x) - cx;
                real_t dy = static_cast<real_t>(y) - cy;
                real_t dz = static_cast<real_t>(z) - cz;

                if (dx * dx * inv_rx2 + dy * dy * inv_ry2 + dz * dz * inv_rz2 <= 1.0f) {
                    grid.set_cell_type(x, y, z, CellType::Solid);
                }
            }
        }
    }
}

} // namespace wake
