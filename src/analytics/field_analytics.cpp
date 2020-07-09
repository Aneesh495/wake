#include "wake/analytics/field_analytics.hpp"
#include <cmath>
#include <algorithm>

namespace ispc {
extern "C" {
    void compute_vorticity_2d(
        const float ux[], const float uy[], float vorticity[], int nx, int ny);
    void compute_vorticity_and_q_criterion_3d(
        const float ux[], const float uy[], const float uz[],
        float omega_x[], float omega_y[], float omega_z[],
        float q_criterion[], int nx, int ny, int nz);
}
}

namespace wake {

void FieldAnalytics::compute_vorticity(const LatticeGrid2D& grid, AlignedBuffer<real_t>& vorticity_out) {
    if (vorticity_out.size() != grid.total_nodes()) {
        vorticity_out.resize(grid.total_nodes());
    }
    ispc::compute_vorticity_2d(
        grid.ux_data(), grid.uy_data(), vorticity_out.data(), grid.nx(), grid.ny()
    );
}

void FieldAnalytics::compute_3d_features(
    const LatticeGrid3D& grid,
    AlignedBuffer<real_t>& omega_x_out,
    AlignedBuffer<real_t>& omega_y_out,
    AlignedBuffer<real_t>& omega_z_out,
    AlignedBuffer<real_t>& q_criterion_out
) {
    std::size_t n = grid.total_nodes();
    if (omega_x_out.size() != n) omega_x_out.resize(n);
    if (omega_y_out.size() != n) omega_y_out.resize(n);
    if (omega_z_out.size() != n) omega_z_out.resize(n);
    if (q_criterion_out.size() != n) q_criterion_out.resize(n);

    ispc::compute_vorticity_and_q_criterion_3d(
        grid.ux_data(), grid.uy_data(), grid.uz_data(),
        omega_x_out.data(), omega_y_out.data(), omega_z_out.data(),
        q_criterion_out.data(),
        grid.nx(), grid.ny(), grid.nz()
    );
}

real_t FieldAnalytics::compute_max_vorticity_2d(const AlignedBuffer<real_t>& vorticity) {
    real_t max_v = 0.0f;
    for (std::size_t i = 0; i < vorticity.size(); ++i) {
        max_v = std::max(max_v, std::abs(vorticity[i]));
    }
    return max_v;
}

} // namespace wake
