#pragma once

#include "wake/core/types.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include "wake/export/terminal_renderer.hpp"
#include <string>
#include <vector>

namespace wake {

class ImageExporter {
public:
    // Export 2D grid field as binary PPM (P6) image
    static bool export_ppm_2d(
        const std::string& filepath,
        const LatticeGrid2D& grid,
        const real_t* scalar_field = nullptr, // if null, uses velocity speed
        Colormap cmap = Colormap::Viridis,
        real_t min_val = 0.0f,
        real_t max_val = 0.15f
    );

    // Export 2D vorticity field as binary PPM (P6) with diverging CoolWarm colormap
    static bool export_ppm_vorticity_2d(
        const std::string& filepath,
        const LatticeGrid2D& grid,
        const real_t* vorticity,
        real_t max_vort = 0.05f
    );

    // Export 3D orthogonal slices (XY, XZ, YZ) to PPM images
    static bool export_ppm_slice_xy(
        const std::string& filepath,
        const LatticeGrid3D& grid,
        int z_slice = -1,
        Colormap cmap = Colormap::Viridis,
        real_t min_val = 0.0f,
        real_t max_val = 0.15f
    );

    static bool export_ppm_slice_xz(
        const std::string& filepath,
        const LatticeGrid3D& grid,
        int y_slice = -1,
        Colormap cmap = Colormap::Viridis,
        real_t min_val = 0.0f,
        real_t max_val = 0.15f
    );

private:
    static RGBColor map_color(Colormap cmap, real_t norm);
};

} // namespace wake
