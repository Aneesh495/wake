#pragma once

#include "wake/core/types.hpp"
#include "wake/core/memory.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/grid/lattice_grid_3d.hpp"

namespace wake {

// High level flow field analytics engine
class FieldAnalytics {
public:
    // Compute 2D vorticity field for a given 2D lattice grid
    static void compute_vorticity(const LatticeGrid2D& grid, AlignedBuffer<real_t>& vorticity_out);

    // Compute 3D vorticity vector and Q-criterion for a given 3D lattice grid
    static void compute_3d_features(
        const LatticeGrid3D& grid,
        AlignedBuffer<real_t>& omega_x_out,
        AlignedBuffer<real_t>& omega_y_out,
        AlignedBuffer<real_t>& omega_z_out,
        AlignedBuffer<real_t>& q_criterion_out
    );

    // Compute maximum vorticity magnitude in domain
    static real_t compute_max_vorticity_2d(const AlignedBuffer<real_t>& vorticity);
};

} // namespace wake
