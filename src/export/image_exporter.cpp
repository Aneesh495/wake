#include "wake/export/image_exporter.hpp"
#include <fstream>
#include <algorithm>
#include <cmath>

namespace wake {

RGBColor ImageExporter::map_color(Colormap cmap, real_t norm) {
    TerminalRenderer temp_renderer(10, 10, cmap);
    return temp_renderer.map_value(norm);
}

bool ImageExporter::export_ppm_2d(
    const std::string& filepath,
    const LatticeGrid2D& grid,
    const real_t* scalar_field,
    Colormap cmap,
    real_t min_val,
    real_t max_val
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }

    int nx = grid.nx();
    int ny = grid.ny();

    // P6 header
    out << "P6\n" << nx << " " << ny << "\n255\n";

    real_t val_range = (max_val - min_val > 1e-8f) ? (max_val - min_val) : 1.0f;
    std::vector<uint8_t> buffer(nx * ny * 3);

    const RGBColor solid_color{80, 80, 85};

    for (int y = 0; y < ny; ++y) {
        // Invert Y so standard image orientation has y=0 at bottom
        int grid_y = (ny - 1) - y;
        for (int x = 0; x < nx; ++x) {
            int buf_idx = (y * nx + x) * 3;

            if (grid.get_cell_type(x, grid_y) == CellType::Solid) {
                buffer[buf_idx]     = solid_color.r;
                buffer[buf_idx + 1] = solid_color.g;
                buffer[buf_idx + 2] = solid_color.b;
            } else {
                real_t val = scalar_field ? scalar_field[grid_y * nx + x] : grid.speed(x, grid_y);
                real_t norm = std::clamp((val - min_val) / val_range, 0.0f, 1.0f);
                RGBColor c = map_color(cmap, norm);

                buffer[buf_idx]     = c.r;
                buffer[buf_idx + 1] = c.g;
                buffer[buf_idx + 2] = c.b;
            }
        }
    }

    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    return true;
}

bool ImageExporter::export_ppm_vorticity_2d(
    const std::string& filepath,
    const LatticeGrid2D& grid,
    const real_t* vorticity,
    real_t max_vort
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }

    int nx = grid.nx();
    int ny = grid.ny();

    out << "P6\n" << nx << " " << ny << "\n255\n";

    real_t denom = (max_vort > 1e-8f) ? (2.0f * max_vort) : 1.0f;
    std::vector<uint8_t> buffer(nx * ny * 3);
    const RGBColor solid_color{80, 80, 85};

    for (int y = 0; y < ny; ++y) {
        int grid_y = (ny - 1) - y;
        for (int x = 0; x < nx; ++x) {
            int buf_idx = (y * nx + x) * 3;

            if (grid.get_cell_type(x, grid_y) == CellType::Solid) {
                buffer[buf_idx]     = solid_color.r;
                buffer[buf_idx + 1] = solid_color.g;
                buffer[buf_idx + 2] = solid_color.b;
            } else {
                real_t omega = vorticity ? vorticity[grid_y * nx + x] : 0.0f;
                real_t norm = std::clamp((omega + max_vort) / denom, 0.0f, 1.0f);
                RGBColor c = map_color(Colormap::CoolWarm, norm);

                buffer[buf_idx]     = c.r;
                buffer[buf_idx + 1] = c.g;
                buffer[buf_idx + 2] = c.b;
            }
        }
    }

    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    return true;
}

bool ImageExporter::export_ppm_slice_xy(
    const std::string& filepath,
    const LatticeGrid3D& grid,
    int z_slice,
    Colormap cmap,
    real_t min_val,
    real_t max_val
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }

    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();

    if (z_slice < 0 || z_slice >= nz) {
        z_slice = nz / 2;
    }

    out << "P6\n" << nx << " " << ny << "\n255\n";

    real_t val_range = (max_val - min_val > 1e-8f) ? (max_val - min_val) : 1.0f;
    std::vector<uint8_t> buffer(nx * ny * 3);
    const RGBColor solid_color{80, 80, 85};

    for (int y = 0; y < ny; ++y) {
        int grid_y = (ny - 1) - y;
        for (int x = 0; x < nx; ++x) {
            int buf_idx = (y * nx + x) * 3;

            if (grid.get_cell_type(x, grid_y, z_slice) == CellType::Solid) {
                buffer[buf_idx]     = solid_color.r;
                buffer[buf_idx + 1] = solid_color.g;
                buffer[buf_idx + 2] = solid_color.b;
            } else {
                real_t val = grid.speed(x, grid_y, z_slice);
                real_t norm = std::clamp((val - min_val) / val_range, 0.0f, 1.0f);
                RGBColor c = map_color(cmap, norm);

                buffer[buf_idx]     = c.r;
                buffer[buf_idx + 1] = c.g;
                buffer[buf_idx + 2] = c.b;
            }
        }
    }

    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    return true;
}

bool ImageExporter::export_ppm_slice_xz(
    const std::string& filepath,
    const LatticeGrid3D& grid,
    int y_slice,
    Colormap cmap,
    real_t min_val,
    real_t max_val
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open()) {
        return false;
    }

    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();

    if (y_slice < 0 || y_slice >= ny) {
        y_slice = ny / 2;
    }

    out << "P6\n" << nx << " " << nz << "\n255\n";

    real_t val_range = (max_val - min_val > 1e-8f) ? (max_val - min_val) : 1.0f;
    std::vector<uint8_t> buffer(nx * nz * 3);
    const RGBColor solid_color{80, 80, 85};

    for (int z = 0; z < nz; ++z) {
        int grid_z = (nz - 1) - z;
        for (int x = 0; x < nx; ++x) {
            int buf_idx = (z * nx + x) * 3;

            if (grid.get_cell_type(x, y_slice, grid_z) == CellType::Solid) {
                buffer[buf_idx]     = solid_color.r;
                buffer[buf_idx + 1] = solid_color.g;
                buffer[buf_idx + 2] = solid_color.b;
            } else {
                real_t val = grid.speed(x, y_slice, grid_z);
                real_t norm = std::clamp((val - min_val) / val_range, 0.0f, 1.0f);
                RGBColor c = map_color(cmap, norm);

                buffer[buf_idx]     = c.r;
                buffer[buf_idx + 1] = c.g;
                buffer[buf_idx + 2] = c.b;
            }
        }
    }

    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
    return true;
}

} // namespace wake
