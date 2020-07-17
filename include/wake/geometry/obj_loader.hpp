#pragma once

#include "wake/core/types.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include <string>
#include <vector>

namespace wake {

struct Triangle {
    Vec3 v0;
    Vec3 v1;
    Vec3 v2;
};

class ObjLoader {
public:
    ObjLoader() = default;

    // Load Wavefront OBJ file
    bool load_file(const std::string& filepath);

    // Load from OBJ format string in memory
    bool load_from_string(const std::string& obj_data);

    // Get parsed mesh statistics
    size_t vertex_count() const { return vertices_.size(); }
    size_t triangle_count() const { return triangles_.size(); }

    const std::vector<Vec3>& vertices() const { return vertices_; }
    const std::vector<Triangle>& triangles() const { return triangles_; }

    // Mesh transformations
    void translate(const Vec3& offset);
    void scale(real_t s);
    void scale(const Vec3& s);
    void rotate_x(real_t angle_deg);
    void rotate_y(real_t angle_deg);
    void rotate_z(real_t angle_deg);

    // Compute Axis-Aligned Bounding Box (AABB)
    void get_bounds(Vec3& min_out, Vec3& max_out) const;

    // Normalize and fit mesh inside target bounding box
    void fit_into_box(const Vec3& target_min, const Vec3& target_max);

    // Voxelize mesh into LatticeGrid3D
    // If fill_interior is true, marks interior nodes as solid
    void voxelize_into_grid(LatticeGrid3D& grid, bool fill_interior = true) const;

private:
    std::vector<Vec3> vertices_;
    std::vector<Triangle> triangles_;

    static bool test_triangle_box_overlap(const Triangle& tri, const Vec3& box_center, const Vec3& box_half_size);
};

} // namespace wake
