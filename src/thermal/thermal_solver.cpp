#include "wake/thermal/thermal_solver.hpp"
#include <numeric>
#include <cmath>
#include <algorithm>

namespace ispc {
extern "C" {
    void thermal_collide_d2q9(
        const float g_in[],
        float g_out[],
        const uint8_t cell_type[],
        const float ux[],
        const float uy[],
        float temperature[],
        float tau_T,
        int count
    );

    void thermal_compute_buoyancy_force_2d(
        const float temperature[],
        float gy_body[],
        float t_ref,
        float beta,
        float g_gravity,
        int count
    );

    void thermal_apply_isothermal_wall(
        float g[],
        int nx,
        int ny,
        int y_wall,
        float t_wall
    );

    void d2q9_stream_with_boundaries(
        const float f_post_collision[],
        float f_streamed[],
        const uint8_t cell_type[],
        int nx,
        int ny
    );
}
}

namespace wake {

ThermalSolver2D::ThermalSolver2D(int nx, int ny, real_t tau_thermal, real_t t_ref, real_t beta)
    : nx_(nx),
      ny_(ny),
      total_nodes_(static_cast<std::size_t>(nx) * ny),
      tau_T_(tau_thermal),
      t_ref_(t_ref),
      beta_(beta),
      g_src_(total_nodes_ * D2Q9::Q),
      g_dst_(total_nodes_ * D2Q9::Q),
      temperature_(total_nodes_),
      gy_buoyancy_(total_nodes_) {
    initialize(t_ref_);
}

void ThermalSolver2D::initialize(real_t t_init) {
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        temperature_[i] = t_init;
        gy_buoyancy_[i] = 0.0f;
        for (int q = 0; q < D2Q9::Q; ++q) {
            g_src_[q * total_nodes_ + i] = D2Q9::W[q] * t_init;
            g_dst_[q * total_nodes_ + i] = D2Q9::W[q] * t_init;
        }
    }
}

void ThermalSolver2D::set_isothermal_wall_south(real_t t_hot) {
    has_hot_south_ = true;
    t_hot_ = t_hot;
}

void ThermalSolver2D::set_isothermal_wall_north(real_t t_cold) {
    has_cold_north_ = true;
    t_cold_ = t_cold;
}

void ThermalSolver2D::step(LatticeGrid2D& fluid_grid, real_t g_gravity) {
    // 1. Collide thermal distributions
    ispc::thermal_collide_d2q9(
        g_src_.data(),
        g_dst_.data(),
        fluid_grid.cell_type_data(),
        fluid_grid.ux_data(),
        fluid_grid.uy_data(),
        temperature_.data(),
        tau_T_,
        static_cast<int>(total_nodes_)
    );

    // 2. Isothermal wall boundary conditions on post-collision distributions
    if (has_hot_south_) {
        ispc::thermal_apply_isothermal_wall(g_dst_.data(), nx_, ny_, 0, t_hot_);
    }
    if (has_cold_north_) {
        ispc::thermal_apply_isothermal_wall(g_dst_.data(), nx_, ny_, ny_ - 1, t_cold_);
    }

    // 3. Pull streaming for thermal scalar distributions
    ispc::d2q9_stream_with_boundaries(
        g_dst_.data(),
        g_src_.data(),
        fluid_grid.cell_type_data(),
        nx_,
        ny_
    );

    // 4. Compute Boussinesq thermal buoyancy acceleration
    std::fill(gy_buoyancy_.data(), gy_buoyancy_.data() + total_nodes_, 0.0f);
    ispc::thermal_compute_buoyancy_force_2d(
        temperature_.data(),
        gy_buoyancy_.data(),
        t_ref_,
        beta_,
        g_gravity,
        static_cast<int>(total_nodes_)
    );
}

real_t ThermalSolver2D::compute_mean_temperature() const {
    double sum = 0.0;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        sum += temperature_[i];
    }
    return static_cast<real_t>(sum / total_nodes_);
}

real_t ThermalSolver2D::compute_nusselt_number() const {
    if (std::abs(t_hot_ - t_cold_) < 1e-6f || ny_ <= 1) {
        return 1.0f;
    }
    // Wall thermal gradient at bottom wall
    double grad_t_sum = 0.0;
    for (int x = 0; x < nx_; ++x) {
        real_t t0 = temperature_[0 * nx_ + x];
        real_t t1 = temperature_[1 * nx_ + x];
        grad_t_sum += (t0 - t1);
    }
    double mean_grad = grad_t_sum / nx_;
    double pure_conduction_grad = (t_hot_ - t_cold_) / (ny_ - 1);
    return static_cast<real_t>(mean_grad / pure_conduction_grad);
}

} // namespace wake
