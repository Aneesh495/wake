#include "wake/geometry/obj_loader.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <queue>

namespace wake {

bool ObjLoader::load_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        return false;
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return load_from_string(buffer.str());
}

bool ObjLoader::load_from_string(const std::string& obj_data) {
    vertices_.clear();
    triangles_.clear();

    std::stringstream ss(obj_data);
    std::string line;

    while (std::getline(ss, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::stringstream line_stream(line);
        std::string prefix;
        line_stream >> prefix;

        if (prefix == "v") {
            real_t x, y, z;
            if (line_stream >> x >> y >> z) {
                vertices_.push_back(Vec3(x, y, z));
            }
        } else if (prefix == "f") {
            std::vector<int> face_indices;
            std::string vertex_token;
            while (line_stream >> vertex_token) {
                // Parse v, v/vt, v/vt/vn, or v//vn
                size_t slash_pos = vertex_token.find('/');
                int idx = 0;
                if (slash_pos == std::string::npos) {
                    idx = std::stoi(vertex_token);
                } else {
                    idx = std::stoi(vertex_token.substr(0, slash_pos));
                }

                // OBJ indices are 1-based, or negative relative
                if (idx > 0) {
                    idx -= 1;
                } else if (idx < 0) {
                    idx = static_cast<int>(vertices_.size()) + idx;
                }
                face_indices.push_back(idx);
            }

            // Triangulate polygon fans (f 0 1 2, f 0 2 3, ...)
            if (face_indices.size() >= 3) {
                for (size_t i = 1; i + 1 < face_indices.size(); ++i) {
                    int i0 = face_indices[0];
                    int i1 = face_indices[i];
                    int i2 = face_indices[i + 1];

                    if (i0 >= 0 && i0 < static_cast<int>(vertices_.size()) &&
                        i1 >= 0 && i1 < static_cast<int>(vertices_.size()) &&
                        i2 >= 0 && i2 < static_cast<int>(vertices_.size())) {
                        triangles_.push_back({vertices_[i0], vertices_[i1], vertices_[i2]});
                    }
                }
            }
        }
    }

    return !triangles_.empty();
}

void ObjLoader::translate(const Vec3& offset) {
    for (auto& v : vertices_) {
        v = v + offset;
    }
    for (auto& tri : triangles_) {
        tri.v0 = tri.v0 + offset;
        tri.v1 = tri.v1 + offset;
        tri.v2 = tri.v2 + offset;
    }
}

void ObjLoader::scale(real_t s) {
    scale(Vec3(s, s, s));
}

void ObjLoader::scale(const Vec3& s) {
    for (auto& v : vertices_) {
        v.x *= s.x;
        v.y *= s.y;
        v.z *= s.z;
    }
    for (auto& tri : triangles_) {
        tri.v0.x *= s.x; tri.v0.y *= s.y; tri.v0.z *= s.z;
        tri.v1.x *= s.x; tri.v1.y *= s.y; tri.v1.z *= s.z;
        tri.v2.x *= s.x; tri.v2.y *= s.y; tri.v2.z *= s.z;
    }
}

void ObjLoader::rotate_x(real_t angle_deg) {
    real_t rad = angle_deg * (PI / 180.0f);
    real_t cos_a = std::cos(rad);
    real_t sin_a = std::sin(rad);

    auto rot = [cos_a, sin_a](const Vec3& p) {
        return Vec3(p.x, p.y * cos_a - p.z * sin_a, p.y * sin_a + p.z * cos_a);
    };

    for (auto& v : vertices_) v = rot(v);
    for (auto& tri : triangles_) {
        tri.v0 = rot(tri.v0);
        tri.v1 = rot(tri.v1);
        tri.v2 = rot(tri.v2);
    }
}

void ObjLoader::rotate_y(real_t angle_deg) {
    real_t rad = angle_deg * (PI / 180.0f);
    real_t cos_a = std::cos(rad);
    real_t sin_a = std::sin(rad);

    auto rot = [cos_a, sin_a](const Vec3& p) {
        return Vec3(p.x * cos_a + p.z * sin_a, p.y, -p.x * sin_a + p.z * cos_a);
    };

    for (auto& v : vertices_) v = rot(v);
    for (auto& tri : triangles_) {
        tri.v0 = rot(tri.v0);
        tri.v1 = rot(tri.v1);
        tri.v2 = rot(tri.v2);
    }
}

void ObjLoader::rotate_z(real_t angle_deg) {
    real_t rad = angle_deg * (PI / 180.0f);
    real_t cos_a = std::cos(rad);
    real_t sin_a = std::sin(rad);

    auto rot = [cos_a, sin_a](const Vec3& p) {
        return Vec3(p.x * cos_a - p.y * sin_a, p.x * sin_a + p.y * cos_a, p.z);
    };

    for (auto& v : vertices_) v = rot(v);
    for (auto& tri : triangles_) {
        tri.v0 = rot(tri.v0);
        tri.v1 = rot(tri.v1);
        tri.v2 = rot(tri.v2);
    }
}

void ObjLoader::get_bounds(Vec3& min_out, Vec3& max_out) const {
    if (vertices_.empty()) {
        min_out = Vec3(0, 0, 0);
        max_out = Vec3(0, 0, 0);
        return;
    }

    min_out = vertices_[0];
    max_out = vertices_[0];

    for (const auto& v : vertices_) {
        min_out.x = std::min(min_out.x, v.x);
        min_out.y = std::min(min_out.y, v.y);
        min_out.z = std::min(min_out.z, v.z);

        max_out.x = std::max(max_out.x, v.x);
        max_out.y = std::max(max_out.y, v.y);
        max_out.z = std::max(max_out.z, v.z);
    }
}

void ObjLoader::fit_into_box(const Vec3& target_min, const Vec3& target_max) {
    Vec3 cur_min, cur_max;
    get_bounds(cur_min, cur_max);

    Vec3 cur_size = cur_max - cur_min;
    Vec3 target_size = target_max - target_min;

    if (cur_size.x < 1e-8f || cur_size.y < 1e-8f || cur_size.z < 1e-8f) {
        return;
    }

    real_t sx = target_size.x / cur_size.x;
    real_t sy = target_size.y / cur_size.y;
    real_t sz = target_size.z / cur_size.z;
    real_t uniform_scale = std::min({sx, sy, sz});

    translate(cur_min * -1.0f);
    scale(uniform_scale);

    Vec3 center_offset = (target_size - (cur_size * uniform_scale)) * 0.5f;
    translate(target_min + center_offset);
}

// Separating Axis Theorem (Akenine-Möller 2001)
bool ObjLoader::test_triangle_box_overlap(const Triangle& tri, const Vec3& box_center, const Vec3& box_half_size) {
    // Translate triangle so box center is at origin
    Vec3 v0 = tri.v0 - box_center;
    Vec3 v1 = tri.v1 - box_center;
    Vec3 v2 = tri.v2 - box_center;

    // Triangle edges
    Vec3 e0 = v1 - v0;
    Vec3 e1 = v2 - v1;
    Vec3 e2 = v0 - v2;

    // Test 1: Box normal axes (AABB bounds check)
    real_t min_x = std::min({v0.x, v1.x, v2.x});
    real_t max_x = std::max({v0.x, v1.x, v2.x});
    if (min_x > box_half_size.x || max_x < -box_half_size.x) return false;

    real_t min_y = std::min({v0.y, v1.y, v2.y});
    real_t max_y = std::max({v0.y, v1.y, v2.y});
    if (min_y > box_half_size.y || max_y < -box_half_size.y) return false;

    real_t min_z = std::min({v0.z, v1.z, v2.z});
    real_t max_z = std::max({v0.z, v1.z, v2.z});
    if (min_z > box_half_size.z || max_z < -box_half_size.z) return false;

    // Test 2: Triangle plane normal
    Vec3 normal = e0.cross(e1);
    real_t d = normal.dot(v0);

    Vec3 vmin, vmax;
    if (normal.x > 0.0f) { vmin.x = -box_half_size.x; vmax.x = box_half_size.x; }
    else { vmin.x = box_half_size.x; vmax.x = -box_half_size.x; }

    if (normal.y > 0.0f) { vmin.y = -box_half_size.y; vmax.y = box_half_size.y; }
    else { vmin.y = box_half_size.y; vmax.y = -box_half_size.y; }

    if (normal.z > 0.0f) { vmin.z = -box_half_size.z; vmax.z = box_half_size.z; }
    else { vmin.z = box_half_size.z; vmax.z = -box_half_size.z; }

    if (normal.dot(vmin) > d || normal.dot(vmax) < d) return false;

    // Test 3: 9 cross-product axes
    auto test_axis = [&](real_t a, real_t b, real_t fa, real_t fb, real_t v0_1, real_t v0_2, real_t v1_1, real_t v1_2, real_t v2_1, real_t v2_2, real_t h1, real_t h2) {
        real_t p0 = a * v0_1 - b * v0_2;
        real_t p1 = a * v1_1 - b * v1_2;
        real_t p2 = a * v2_1 - b * v2_2;
        real_t min_p = std::min({p0, p1, p2});
        real_t max_p = std::max({p0, p1, p2});
        real_t rad = fa * h1 + fb * h2;
        return (min_p <= rad && max_p >= -rad);
    };

    real_t fe0x = std::abs(e0.x);
    real_t fe0y = std::abs(e0.y);
    real_t fe0z = std::abs(e0.z);

    if (!test_axis(e0.z, e0.y, fe0z, fe0y, v0.y, v0.z, v1.y, v1.z, v2.y, v2.z, box_half_size.y, box_half_size.z)) return false;
    if (!test_axis(-e0.z, e0.x, fe0z, fe0x, v0.x, v0.z, v1.x, v1.z, v2.x, v2.z, box_half_size.x, box_half_size.z)) return false;
    if (!test_axis(e0.y, -e0.x, fe0y, fe0x, v0.x, v0.y, v1.x, v1.y, v2.x, v2.y, box_half_size.x, box_half_size.y)) return false;

    real_t fe1x = std::abs(e1.x);
    real_t fe1y = std::abs(e1.y);
    real_t fe1z = std::abs(e1.z);

    if (!test_axis(e1.z, e1.y, fe1z, fe1y, v0.y, v0.z, v1.y, v1.z, v2.y, v2.z, box_half_size.y, box_half_size.z)) return false;
    if (!test_axis(-e1.z, e1.x, fe1z, fe1x, v0.x, v0.z, v1.x, v1.z, v2.x, v2.z, box_half_size.x, box_half_size.z)) return false;
    if (!test_axis(e1.y, -e1.x, fe1y, fe1x, v0.x, v0.y, v1.x, v1.y, v2.x, v2.y, box_half_size.x, box_half_size.y)) return false;

    real_t fe2x = std::abs(e2.x);
    real_t fe2y = std::abs(e2.y);
    real_t fe2z = std::abs(e2.z);

    if (!test_axis(e2.z, e2.y, fe2z, fe2y, v0.y, v0.z, v1.y, v1.z, v2.y, v2.z, box_half_size.y, box_half_size.z)) return false;
    if (!test_axis(-e2.z, e2.x, fe2z, fe2x, v0.x, v0.z, v1.x, v1.z, v2.x, v2.z, box_half_size.x, box_half_size.z)) return false;
    if (!test_axis(e2.y, -e2.x, fe2y, fe2x, v0.x, v0.y, v1.x, v1.y, v2.x, v2.y, box_half_size.x, box_half_size.y)) return false;

    return true;
}

void ObjLoader::voxelize_into_grid(LatticeGrid3D& grid, bool fill_interior) const {
    const Vec3 box_half(0.5f, 0.5f, 0.5f);

    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();

    // Step 1: Surface voxelization using SAT
    for (const auto& tri : triangles_) {
        real_t min_x = std::min({tri.v0.x, tri.v1.x, tri.v2.x});
        real_t max_x = std::max({tri.v0.x, tri.v1.x, tri.v2.x});
        real_t min_y = std::min({tri.v0.y, tri.v1.y, tri.v2.y});
        real_t max_y = std::max({tri.v0.y, tri.v1.y, tri.v2.y});
        real_t min_z = std::min({tri.v0.z, tri.v1.z, tri.v2.z});
        real_t max_z = std::max({tri.v0.z, tri.v1.z, tri.v2.z});

        int i_min = std::max(0, static_cast<int>(std::floor(min_x)));
        int i_max = std::min(nx - 1, static_cast<int>(std::ceil(max_x)));
        int j_min = std::max(0, static_cast<int>(std::floor(min_y)));
        int j_max = std::min(ny - 1, static_cast<int>(std::ceil(max_y)));
        int k_min = std::max(0, static_cast<int>(std::floor(min_z)));
        int k_max = std::min(nz - 1, static_cast<int>(std::ceil(max_z)));

        for (int k = k_min; k <= k_max; ++k) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int i = i_min; i <= i_max; ++i) {
                    Vec3 center(static_cast<real_t>(i) + 0.5f,
                                static_cast<real_t>(j) + 0.5f,
                                static_cast<real_t>(k) + 0.5f);

                    if (test_triangle_box_overlap(tri, center, box_half)) {
                        grid.set_cell_type(i, j, k, CellType::Solid);
                    }
                }
            }
        }
    }

    if (!fill_interior) {
        return;
    }

    // Step 2: Flood fill from exterior boundary to find all exterior fluid nodes
    std::vector<uint8_t> visited(nx * ny * nz, 0);
    std::queue<int> q;

    auto push_node = [&](int x, int y, int z) {
        int idx = (z * ny + y) * nx + x;
        if (!visited[idx] && grid.get_cell_type(x, y, z) != CellType::Solid) {
            visited[idx] = 1;
            q.push(idx);
        }
    };

    // Seed outer boundary faces
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            push_node(0, j, k);
            push_node(nx - 1, j, k);
        }
    }
    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
            push_node(i, 0, k);
            push_node(i, ny - 1, k);
        }
    }
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            push_node(i, j, 0);
            push_node(i, j, nz - 1);
        }
    }

    const int dx[] = {1, -1, 0, 0, 0, 0};
    const int dy[] = {0, 0, 1, -1, 0, 0};
    const int dz[] = {0, 0, 0, 0, 1, -1};

    while (!q.empty()) {
        int curr = q.front();
        q.pop();

        int x = curr % nx;
        int y = (curr / nx) % ny;
        int z = curr / (nx * ny);

        for (int d = 0; d < 6; ++d) {
            int nx_pos = x + dx[d];
            int ny_pos = y + dy[d];
            int nz_pos = z + dz[d];

            if (nx_pos >= 0 && nx_pos < nx &&
                ny_pos >= 0 && ny_pos < ny &&
                nz_pos >= 0 && nz_pos < nz) {
                int n_idx = (nz_pos * ny + ny_pos) * nx + nx_pos;
                if (!visited[n_idx] && grid.get_cell_type(nx_pos, ny_pos, nz_pos) != CellType::Solid) {
                    visited[n_idx] = 1;
                    q.push(n_idx);
                }
            }
        }
    }

    // Any cell not reachable from outside is interior solid
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = (k * ny + j) * nx + i;
                if (!visited[idx]) {
                    grid.set_cell_type(i, j, k, CellType::Solid);
                }
            }
        }
    }
}

} // namespace wake
