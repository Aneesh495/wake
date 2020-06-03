#include "wake/grid/lattice_grid_3d.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

// C linkage declarations for ISPC 3D functions
namespace ispc {
extern "C" {
    void d3q19_compute_macroscopic_moments(
        const float f[], float rho[], float ux[], float uy[], float uz[], int count);
    void d3q19_compute_macroscopic_moments_forced(
        const float f[], float rho[], float ux[], float uy[], float uz[],
        float gx, float gy, float gz, int count);
    void d3q19_compute_velocity_magnitude(
        const float ux[], const float uy[], const float uz[], float speed[], int count);
    void d3q19_collide_bgk(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float uz[],
        float tau, float gx, float gy, float gz, int count);
    void d3q19_collide_trt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float uz[],
        float tau_s, float tau_a, float gx, float gy, float gz, int count);
    void d3q19_collide_mrt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float uz[],
        float tau_shear, float tau_bulk, float gx, float gy, float gz, int count);
    void d3q19_stream_periodic(
        const float f_src[], float f_dst[], int nx, int ny, int nz);
    void d3q19_stream_with_boundaries(
        const float f_src[], float f_dst[], const uint8_t cell_type[],
        int nx, int ny, int nz, float uw_x, float uw_y, float uw_z);
    void d3q19_compute_smagorinsky_tau(
        const float f[], const float rho[], const float ux[], const float uy[], const float uz[],
        float tau_eff[], float tau_0, float smagorinsky_constant, int count);
    void d3q19_apply_curved_bfl(
        float f_dst[], const float f_src[], const int fluid_indices[],
        const int q_dirs[], const float delta_fractions[],
        const float uw_x[], const float uw_y[], const float uw_z[],
        int nx, int ny, int nz, int num_boundary_links);
}
}

namespace wake {

LatticeGrid3D::LatticeGrid3D(int nx, int ny, int nz, real_t relaxation_time)
    : nx_(nx), ny_(ny), nz_(nz),
      total_nodes_(static_cast<std::size_t>(nx) * ny * nz),
      tau_(relaxation_time) {
    const std::size_t total_q_nodes = D3Q19::Q * total_nodes_;

    f_src_.resize(total_q_nodes);
    f_dst_.resize(total_q_nodes);

    rho_.resize(total_nodes_);
    ux_.resize(total_nodes_);
    uy_.resize(total_nodes_);
    uz_.resize(total_nodes_);
    speed_.resize(total_nodes_);
    tau_eff_.resize(total_nodes_);
    cell_type_.resize(total_nodes_);

    initialize(RHO_0, 0.0f, 0.0f, 0.0f);
}

void LatticeGrid3D::initialize(real_t rho_init, real_t ux_init, real_t uy_init, real_t uz_init) {
    step_count_ = 0;
    real_t u_sq = ux_init * ux_init + uy_init * uy_init + uz_init * uz_init;

    for (int q = 0; q < D3Q19::Q; ++q) {
        real_t ci_u = static_cast<real_t>(D3Q19::CX[q]) * ux_init +
                     static_cast<real_t>(D3Q19::CY[q]) * uy_init +
                     static_cast<real_t>(D3Q19::CZ[q]) * uz_init;
        real_t feq = D3Q19::W[q] * rho_init * (
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
        uz_[i] = uz_init;
        speed_[i] = std::sqrt(u_sq);
        tau_eff_[i] = tau_;
        cell_type_[i] = static_cast<uint8_t>(CellType::Fluid);
    }
}

void LatticeGrid3D::initialize_custom(
    const std::function<void(int x, int y, int z, real_t& rho, real_t& ux, real_t& uy, real_t& uz)>& init_fn
) {
    step_count_ = 0;
    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                std::size_t idx = (static_cast<std::size_t>(z) * ny_ + y) * nx_ + x;
                real_t local_rho = RHO_0;
                real_t local_ux = 0.0f;
                real_t local_uy = 0.0f;
                real_t local_uz = 0.0f;

                init_fn(x, y, z, local_rho, local_ux, local_uy, local_uz);

                rho_[idx] = local_rho;
                ux_[idx] = local_ux;
                uy_[idx] = local_uy;
                uz_[idx] = local_uz;
                real_t u_sq = local_ux * local_ux + local_uy * local_uy + local_uz * local_uz;
                speed_[idx] = std::sqrt(u_sq);
                tau_eff_[idx] = tau_;

                for (int q = 0; q < D3Q19::Q; ++q) {
                    real_t ci_u = static_cast<real_t>(D3Q19::CX[q]) * local_ux +
                                 static_cast<real_t>(D3Q19::CY[q]) * local_uy +
                                 static_cast<real_t>(D3Q19::CZ[q]) * local_uz;
                    real_t feq = D3Q19::W[q] * local_rho * (
                        1.0f + INV_CS2 * ci_u + INV_2_CS4 * ci_u * ci_u - INV_2_CS2 * u_sq
                    );
                    f_src_[q * total_nodes_ + idx] = feq;
                    f_dst_[q * total_nodes_ + idx] = feq;
                }
            }
        }
    }
}

void LatticeGrid3D::set_cell_type(int x, int y, int z, CellType type) {
    if (x >= 0 && x < nx_ && y >= 0 && y < ny_ && z >= 0 && z < nz_) {
        cell_type_[(static_cast<std::size_t>(z) * ny_ + y) * nx_ + x] = static_cast<uint8_t>(type);
    }
}

CellType LatticeGrid3D::get_cell_type(int x, int y, int z) const {
    if (x >= 0 && x < nx_ && y >= 0 && y < ny_ && z >= 0 && z < nz_) {
        return static_cast<CellType>(cell_type_[(static_cast<std::size_t>(z) * ny_ + y) * nx_ + x]);
    }
    return CellType::Solid;
}

void LatticeGrid3D::set_solid_box(int x0, int y0, int z0, int x1, int y1, int z1) {
    for (int z = std::max(0, z0); z <= std::min(nz_ - 1, z1); ++z) {
        for (int y = std::max(0, y0); y <= std::min(ny_ - 1, y1); ++y) {
            for (int x = std::max(0, x0); x <= std::min(nx_ - 1, x1); ++x) {
                set_cell_type(x, y, z, CellType::Solid);
            }
        }
    }
}

void LatticeGrid3D::set_solid_sphere(real_t cx, real_t cy, real_t cz, real_t radius) {
    real_t r_sq = radius * radius;
    int min_x = std::max(0, static_cast<int>(cx - radius - 1));
    int max_x = std::min(nx_ - 1, static_cast<int>(cx + radius + 1));
    int min_y = std::max(0, static_cast<int>(cy - radius - 1));
    int max_y = std::min(ny_ - 1, static_cast<int>(cy + radius + 1));
    int min_z = std::max(0, static_cast<int>(cz - radius - 1));
    int max_z = std::min(nz_ - 1, static_cast<int>(cz + radius + 1));

    for (int z = min_z; z <= max_z; ++z) {
        for (int y = min_y; y <= max_y; ++y) {
            for (int x = min_x; x <= max_x; ++x) {
                real_t dx = static_cast<real_t>(x) - cx;
                real_t dy = static_cast<real_t>(y) - cy;
                real_t dz = static_cast<real_t>(z) - cz;
                if (dx * dx + dy * dy + dz * dz <= r_sq) {
                    set_cell_type(x, y, z, CellType::Solid);
                }
            }
        }
    }
}

void LatticeGrid3D::set_solid_cylinder_z(real_t cx, real_t cy, real_t radius) {
    real_t r_sq = radius * radius;
    int min_x = std::max(0, static_cast<int>(cx - radius - 1));
    int max_x = std::min(nx_ - 1, static_cast<int>(cx + radius + 1));
    int min_y = std::max(0, static_cast<int>(cy - radius - 1));
    int max_y = std::min(ny_ - 1, static_cast<int>(cy + radius + 1));

    for (int z = 0; z < nz_; ++z) {
        for (int y = min_y; y <= max_y; ++y) {
            for (int x = min_x; x <= max_x; ++x) {
                real_t dx = static_cast<real_t>(x) - cx;
                real_t dy = static_cast<real_t>(y) - cy;
                if (dx * dx + dy * dy <= r_sq) {
                    set_cell_type(x, y, z, CellType::Solid);
                }
            }
        }
    }
}

void LatticeGrid3D::set_moving_wall_top(real_t u_wall_x, real_t u_wall_z) {
    has_moving_wall_top_ = true;
    u_wall_top_x_ = u_wall_x;
    u_wall_top_z_ = u_wall_z;
    for (int z = 0; z < nz_; ++z) {
        for (int x = 0; x < nx_; ++x) {
            set_cell_type(x, ny_ - 1, z, CellType::MovingWall);
        }
    }
}

void LatticeGrid3D::build_curved_sphere_links(real_t cx, real_t cy, real_t cz, real_t radius) {
    std::vector<int> fluid_indices;
    std::vector<int> q_dirs;
    std::vector<real_t> delta_fractions;
    std::vector<real_t> uw_x, uw_y, uw_z;

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                if (get_cell_type(x, y, z) != CellType::Fluid) continue;
                int idx_f = (z * ny_ + y) * nx_ + x;

                for (int q = 1; q < D3Q19::Q; ++q) {
                    int nx_x = (x + D3Q19::CX[q] + nx_) % nx_;
                    int ny_y = y + D3Q19::CY[q];
                    int nz_z = z + D3Q19::CZ[q];

                    if (ny_y >= 0 && ny_y < ny_ && nz_z >= 0 && nz_z < nz_) {
                        if (get_cell_type(nx_x, ny_y, nz_z) == CellType::Solid) {
                            // Find analytical intersection of line segment with sphere
                            Vec3 p0(static_cast<real_t>(x) - cx,
                                    static_cast<real_t>(y) - cy,
                                    static_cast<real_t>(z) - cz);
                            Vec3 d(static_cast<real_t>(D3Q19::CX[q]),
                                   static_cast<real_t>(D3Q19::CY[q]),
                                   static_cast<real_t>(D3Q19::CZ[q]));

                            // Solve |p0 + delta * d|^2 = radius^2
                            real_t a = d.norm_sq();
                            real_t b = 2.0f * p0.dot(d);
                            real_t c = p0.norm_sq() - radius * radius;
                            real_t disc = b * b - 4.0f * a * c;

                            real_t delta = 0.5f;
                            if (disc >= 0.0f) {
                                real_t t1 = (-b - std::sqrt(disc)) / (2.0f * a);
                                real_t t2 = (-b + std::sqrt(disc)) / (2.0f * a);
                                if (t1 >= 0.0f && t1 <= 1.0f) delta = t1;
                                else if (t2 >= 0.0f && t2 <= 1.0f) delta = t2;
                            }
                            delta = std::clamp(delta, 0.01f, 0.99f);

                            fluid_indices.push_back(idx_f);
                            q_dirs.push_back(q);
                            delta_fractions.push_back(delta);
                            uw_x.push_back(0.0f);
                            uw_y.push_back(0.0f);
                            uw_z.push_back(0.0f);
                        }
                    }
                }
            }
        }
    }

    num_bfl_links_ = static_cast<int>(fluid_indices.size());
    if (num_bfl_links_ > 0) {
        use_curved_bfl_ = true;
        bfl_fluid_indices_.resize(num_bfl_links_);
        bfl_q_dirs_.resize(num_bfl_links_);
        bfl_delta_fractions_.resize(num_bfl_links_);
        bfl_uw_x_.resize(num_bfl_links_);
        bfl_uw_y_.resize(num_bfl_links_);
        bfl_uw_z_.resize(num_bfl_links_);

        for (int i = 0; i < num_bfl_links_; ++i) {
            bfl_fluid_indices_[i] = fluid_indices[i];
            bfl_q_dirs_[i] = q_dirs[i];
            bfl_delta_fractions_[i] = delta_fractions[i];
            bfl_uw_x_[i] = uw_x[i];
            bfl_uw_y_[i] = uw_y[i];
            bfl_uw_z_[i] = uw_z[i];
        }
    }
}

void LatticeGrid3D::update_macroscopic() {
    ispc::d3q19_compute_macroscopic_moments(
        f_src_.data(), rho_.data(), ux_.data(), uy_.data(), uz_.data(), static_cast<int>(total_nodes_)
    );
    ispc::d3q19_compute_velocity_magnitude(
        ux_.data(), uy_.data(), uz_.data(), speed_.data(), static_cast<int>(total_nodes_)
    );
}

void LatticeGrid3D::step(
    CollisionType collision,
    real_t gx,
    real_t gy,
    real_t gz,
    bool use_les,
    real_t smag_const
) {
    int count = static_cast<int>(total_nodes_);

    // 1. Optional LES eddy-viscosity update
    real_t effective_tau = tau_;
    if (use_les) {
        ispc::d3q19_compute_smagorinsky_tau(
            f_src_.data(), rho_.data(), ux_.data(), uy_.data(), uz_.data(),
            tau_eff_.data(), tau_, smag_const, count
        );
    }

    // 2. Collision step
    if (collision == CollisionType::BGK) {
        ispc::d3q19_collide_bgk(
            f_src_.data(), f_dst_.data(), cell_type_.data(),
            rho_.data(), ux_.data(), uy_.data(), uz_.data(),
            effective_tau, gx, gy, gz, count
        );
    } else if (collision == CollisionType::TRT) {
        real_t tau_s = effective_tau;
        real_t diff = std::max(1e-4f, tau_s - 0.5f);
        real_t tau_a = 0.5f + 0.25f / diff;
        ispc::d3q19_collide_trt(
            f_src_.data(), f_dst_.data(), cell_type_.data(),
            rho_.data(), ux_.data(), uy_.data(), uz_.data(),
            tau_s, tau_a, gx, gy, gz, count
        );
    } else { // MRT
        real_t tau_shear = effective_tau;
        real_t tau_bulk = 1.19f;
        ispc::d3q19_collide_mrt(
            f_src_.data(), f_dst_.data(), cell_type_.data(),
            rho_.data(), ux_.data(), uy_.data(), uz_.data(),
            tau_shear, tau_bulk, gx, gy, gz, count
        );
    }

    // 3. Streaming step with integrated halfway bounce-back
    ispc::d3q19_stream_with_boundaries(
        f_dst_.data(), f_src_.data(), cell_type_.data(),
        nx_, ny_, nz_,
        has_moving_wall_top_ ? u_wall_top_x_ : 0.0f,
        0.0f,
        has_moving_wall_top_ ? u_wall_top_z_ : 0.0f
    );

    // 4. Second-order Bouzidi curved boundary correction if active
    if (use_curved_bfl_ && num_bfl_links_ > 0) {
        ispc::d3q19_apply_curved_bfl(
            f_src_.data(), f_dst_.data(),
            bfl_fluid_indices_.data(), bfl_q_dirs_.data(),
            bfl_delta_fractions_.data(),
            bfl_uw_x_.data(), bfl_uw_y_.data(), bfl_uw_z_.data(),
            nx_, ny_, nz_, num_bfl_links_
        );
    }

    // 5. Update macroscopic moments with body force correction
    ispc::d3q19_compute_macroscopic_moments_forced(
        f_src_.data(), rho_.data(), ux_.data(), uy_.data(), uz_.data(),
        gx, gy, gz, count
    );
    ispc::d3q19_compute_velocity_magnitude(
        ux_.data(), uy_.data(), uz_.data(), speed_.data(), count
    );

    ++step_count_;
}

real_t LatticeGrid3D::compute_total_mass() const {
    real_t total = 0.0f;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        if (cell_type_[i] == static_cast<uint8_t>(CellType::Fluid)) {
            total += rho_[i];
        }
    }
    return total;
}

real_t LatticeGrid3D::compute_total_kinetic_energy() const {
    real_t total = 0.0f;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        if (cell_type_[i] == static_cast<uint8_t>(CellType::Fluid)) {
            total += 0.5f * rho_[i] * speed_[i] * speed_[i];
        }
    }
    return total;
}

real_t LatticeGrid3D::compute_max_speed() const {
    real_t max_s = 0.0f;
    for (std::size_t i = 0; i < total_nodes_; ++i) {
        if (cell_type_[i] == static_cast<uint8_t>(CellType::Fluid)) {
            max_s = std::max(max_s, speed_[i]);
        }
    }
    return max_s;
}

} // namespace wake
