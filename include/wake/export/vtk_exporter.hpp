#pragma once

#include "wake/core/types.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include <string>

namespace wake {

class VtkExporter {
public:
    // Export 2D LatticeGrid2D to VTK XML ImageData (.vti) format
    static bool export_vti_2d(
        const std::string& filepath,
        const LatticeGrid2D& grid,
        const real_t* vorticity = nullptr
    );

    // Export 3D LatticeGrid3D to VTK XML ImageData (.vti) format
    static bool export_vti_3d(
        const std::string& filepath,
        const LatticeGrid3D& grid,
        const real_t* q_criterion = nullptr,
        const Vec3* vorticity = nullptr
    );

    // Write a ParaView .pvd time series collection manifest
    static bool write_pvd_manifest(
        const std::string& manifest_path,
        const std::string& file_prefix,
        int start_step,
        int end_step,
        int step_interval,
        real_t dt = 1.0f
    );
};

} // namespace wake
