#include "wake/grid/lattice_grid_2d.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

// C linkage declarations for ISPC exported functions
namespace ispc {
extern "C" {
    void d2q9_compute_macroscopic_moments(
        const float f[], float rho[], float ux[], float uy[], int count);
    void d2q9_compute_velocity_magnitude(
        const float ux[], const float uy[], float speed[], int count);
    void d2q9_collide_bgk(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float tau, float gx, float gy, int count);
    void d2q9_collide_trt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float tau_s, float tau_a, float gx, float gy, int count);
    void d2q9_collide_mrt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float tau_shear, float tau_bulk, float gx, float gy, int count);
    void d2q9_stream_periodic(
        const float f_src[], float f_dst[], int nx, int ny);
    void d2q9_stream_interior(
        const float f_src[], float f_dst[], int nx, int ny);
    void d2q9_apply_bounce_back(
        float f[], const uint8_t cell_type[], float uw_x, float uw_y, int count);
    void d2q9_apply_zou_he_inflow_west(
        float f[], float rho[], float ux[], float uy[], float u_in_x, float u_in_y, int nx, int ny);
    void d2q9_apply_zou_he_outflow_east(
        float f[], float rho[], float ux[], float uy[], float rho_out, int nx, int ny);
    void d2q9_apply_convective_outflow_east(
        float f[], const float f_prev[], float u_conv, int nx, int ny);
    void d2q9_compute_smagorinsky_tau(
        const float f[], const float rho[], const float ux[], const float uy[],
        float tau_eff[], float tau_0, float smagorinsky_constant, int count);
}
}

namespace wake {

LatticeGrid2D::LatticeGrid2D(int nx, int ny, real_t relaxation_time)
    : nx_(nx), ny_(ny), total_nodes_(static_cast<std::size_t>(nx) * ny), tau_(relaxation_time) {
    const std::size_t total_q_nodes = D2Q9::Q * total_nodes_;

    f_src_.resize(total_q_nodes);
    f_dst_.resize(total_q_nodes);
    f_prev_outflow_.resize(total_q_nodes);

    rho_.resize(total_nodes_);
    ux_.resize(total_nodes_);
    uy_.resize(total_nodes_);
    speed_.resize(total_nodes_);
    tau_eff_.resize(total_nodes_);
    cell_type_.resize(total_nodes_);

    initialize(RHO_0, 0.0f, 0.0f);
}

void LatticeGrid2D::initialize(real_t rho_init, real_t ux_init, real_t uy_init) {
    step_count_ = 0;
    real_t u_sq = ux_init * ux_init + uy_init * uy_init;

    for (int q = 0; q < D2Q9::Q; ++q) {
        real_t ci_u = static_cast<real_t>(D2Q9::CX[q]) * ux_init + static_cast<real_t>(D2Q9::CY[q]) * uy_init;
        real_t feq = D2Q9::W[q] * rho_init * (
            1.0f + INV_CS2 * ci_u + INV_2_CS4 * ci_u * ci_u - INV_2_CS2 * u_sq
        );

        for (std::size_t i = 0; i < total_nodes_; ++i) {
            f_src_[q * total_nodes_ + i] = feq;
            f_dst_[q * total_nodes_ + i] = feq;
        }
    }

    for (std::size_t i = 0; i < total_nodes_; ++i) {
        rho_[i] = rho_init;
        ux_[i] = ux_init;
        uy_[i] = uy_init;
        speed_[i] = std::sqrt(u_sq);
        tau_eff_[i] = tau_;
        cell_type_[i] = static_cast<uint8_t>(CellType::Fluid);
    }
}

void LatticeGrid2D::initialize_custom(
    const std::function<void(int x, int y, real_t& rho, real_t& ux, real_t& uy)>& init_fn
) {
    step_count_ = 0;
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            std::size_t idx = static_cast<std::size_t>(y) * nx_ + x;
            real_t local_rho = RHO_0;
            real_t local_ux = 0.0f;
            real_t local_uy = 0.0f;

            init_fn(x, y, local_rho, local_ux, local_uy);

            rho_[idx] = local_rho;
            ux_[idx] = local_ux;
            uy_[idx] = local_uy;
            speed_[idx] = std::sqrt(local_ux * local_ux + local_uy * local_uy);
            tau_eff_[idx] = tau_;

            real_t u_sq = local_ux * local_ux + local_uy * local_uy;
            for (int q = 0; q < D2Q9::Q; ++q) {
                real_t ci_u = static_cast<real_t>(D2Q9::CX[q]) * local_ux + static_cast<real_t>(D2Q9::CY[q]) * local_uy;
                real_t feq = D2Q9::W[q] * local_rho * (
                    1.0f + INV_CS2 * ci_u + INV_2_CS4 * ci_u * ci_u - INV_2_CS2 * u_sq
                );
                f_src_[q * total_nodes_ + idx] = feq;
                f_dst_[q * total_nodes_ + idx] = feq;
            }
        }
    }
}

void LatticeGrid2D::set_cell_type(int x, int y, CellType type) {
    if (x >= 0 && x < nx_ && y >= 0 && y < ny_) {
        cell_type_[static_cast<std::size_t>(y) * nx_ + x] = static_cast<uint8_t>(type);
    }
}

CellType LatticeGrid2D::get_cell_type(int x, int y) const {
    if (x >= 0 && x < nx_ && y >= 0 && y < ny_) {
        return static_cast<CellType>(cell_type_[static_cast<std::size_t>(y) * nx_ + x]);
    }
    return CellType::Solid;
}

void LatticeGrid2D::set_solid_box(int x0, int y0, int x1, int y1) {
    for (int y = std::max(0, y0); y <= std::min(ny_ - 1, y1); ++y) {
        for (int x = std::max(0, x0); x <= std::min(nx_ - 1, x1); ++x) {
            set_cell_type(x, y, CellType::Solid);
        }
    }
}

void LatticeGrid2D::set_solid_circle(real_t cx, real_t cy, real_t radius) {
    real_t r_sq = radius * radius;
    int min_x = std::max(0, static_cast<int>(cx - radius - 1));
    int max_x = std::min(nx_ - 1, static_cast<int>(cx + radius + 1));
    int min_y = std::max(0, static_cast<int>(cy - radius - 1));
    int max_y = std::min(ny_ - 1, static_cast<int>(cy + radius + 1));

    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            real_t dx = static_cast<real_t>(x) - cx;
            real_t dy = static_cast<real_t>(y) - cy;
            if (dx * dx + dy * dy <= r_sq) {
                set_cell_type(x, y, CellType::Solid);
            }
        }
    }
}

void LatticeGrid2D::set_inflow_west(real_t u_in_x, real_t u_in_y) {
    has_inflow_west_ = true;
    u_inflow_x_ = u_in_x;
    u_inflow_y_ = u_in_y;
    for (int y = 0; y < ny_; ++y) {
        set_cell_type(0, y, CellType::Inflow);
    }
}

void LatticeGrid2D::set_outflow_east(real_t rho_out) {
    has_outflow_east_ = true;
    rho_outflow_ = rho_out;
    for (int y = 0; y < ny_; ++y) {
        set_cell_type(nx_ - 1, y, CellType::Outflow);
    }
}

void LatticeGrid2D::set_moving_wall_north(real_t u_wall_x) {
    has_moving_wall_north_ = true;
    u_wall_north_x_ = u_wall_x;
    for (int x = 0; x < nx_; ++x) {
        set_cell_type(x, ny_ - 1, CellType::MovingWall);
    }
}

void LatticeGrid2D::update_macroscopic() {
    ispc::d2q9_compute_macroscopic_moments(
        f_src_.data(), rho_.data(), ux_.data(), uy_.data(), static_cast<int>(total_nodes_)
    );
    ispc::d2q9_compute_velocity_magnitude(
        ux_.data(), uy_.data(), speed_.data(), static_cast<int>(total_nodes_)
    );
}

void LatticeGrid2D::step(
    CollisionType collision,
    real_t gx,
    real_t gy,
    bool use_les,
    real_t smag_const
) {
    int count = static_cast<int>(total_nodes_);

    // 1. Optional LES eddy-viscosity update
    real_t effective_tau = tau_;
    if (use_les) {
        ispc::d2q9_compute_smagorinsky_tau(
            f_src_.data(), rho_.data(), ux_.data(), uy_.data(),
            tau_eff_.data(), tau_, smag_const, count
        );
    }

    // 2. Collision step: f_src -> f_dst (temporary storage)
    if (collision == CollisionType::BGK) {
        ispc::d2q9_collide_bgk(
            f_src_.data(), f_dst_.data(), cell_type_.data(),
            rho_.data(), ux_.data(), uy_.data(), effective_tau, gx, gy, count
        );
    } else if (collision == CollisionType::TRT) {
        // Magic parameter Lambda = 1/4 = 0.25 -> tau_a = 0.5 + 0.25 / (tau_s - 0.5)
        real_t tau_s = effective_tau;
        real_t diff = std::max(1e-4f, tau_s - 0.5f);
        real_t tau_a = 0.5f + 0.25f / diff;
        ispc::d2q9_collide_trt(
            f_src_.data(), f_dst_.data(), cell_type_.data(),
            rho_.data(), ux_.data(), uy_.data(), tau_s, tau_a, gx, gy, count
        );
    } else { // MRT
        real_t tau_shear = effective_tau;
        real_t tau_bulk = 1.63f;
        ispc::d2q9_collide_mrt(
            f_src_.data(), f_dst_.data(), cell_type_.data(),
            rho_.data(), ux_.data(), uy_.data(), tau_shear, tau_bulk, gx, gy, count
        );
    }

    // 3. Streaming step: f_dst (post-collision) -> f_src (streamed)
    if (!has_inflow_west_ && !has_outflow_east_) {
        ispc::d2q9_stream_periodic(f_dst_.data(), f_src_.data(), nx_, ny_);
    } else {
        ispc::d2q9_stream_interior(f_dst_.data(), f_src_.data(), nx_, ny_);
    }

    // 4. Boundary conditions on f_src
    // Solid obstacles bounce-back
    ispc::d2q9_apply_bounce_back(
        f_src_.data(), cell_type_.data(),
        has_moving_wall_north_ ? u_wall_north_x_ : 0.0f,
        0.0f, count
    );

    // Inflow at west
    if (has_inflow_west_) {
        ispc::d2q9_apply_zou_he_inflow_west(
            f_src_.data(), rho_.data(), ux_.data(), uy_.data(),
            u_inflow_x_, u_inflow_y_, nx_, ny_
        );
    }

    // Outflow at east
    if (has_outflow_east_) {
        ispc::d2q9_apply_convective_outflow_east(
            f_src_.data(), f_prev_outflow_.data(), u_inflow_x_ > 0.0f ? u_inflow_x_ : 0.05f, nx_, ny_
        );
        std::memcpy(f_prev_outflow_.data(), f_src_.data(), f_src_.bytes());
    }

    // 5. Update macroscopic moments
    update_macroscopic();

    ++step_count_;
}

real_t LatticeGrid2D::compute_total_mass() const {
    real_t total = 0.0f;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        if (cell_type_[i] == static_cast<uint8_t>(CellType::Fluid)) {
            total += rho_[i];
        }
    }
    return total;
}

real_t LatticeGrid2D::compute_total_kinetic_energy() const {
    real_t total = 0.0f;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        if (cell_type_[i] == static_cast<uint8_t>(CellType::Fluid)) {
            total += 0.5f * rho_[i] * speed_[i] * speed_[i];
        }
    }
    return total;
}

real_t LatticeGrid2D::compute_max_speed() const {
    real_t max_s = 0.0f;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        if (cell_type_[i] == static_cast<uint8_t>(CellType::Fluid)) {
            max_s = std::max(max_s, speed_[i]);
        }
    }
    return max_s;
}

} // namespace wake
