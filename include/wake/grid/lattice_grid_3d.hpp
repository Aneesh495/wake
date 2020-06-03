#pragma once

#include "wake/core/types.hpp"
#include "wake/core/constants.hpp"
#include "wake/core/memory.hpp"
#include "wake/core/lattice_models.hpp"
#include <vector>
#include <functional>

namespace wake {

// Curved boundary link descriptor for second-order Bouzidi boundary condition
struct CurvedBoundaryLink {
    int fluid_index;
    int q_dir;
    real_t delta;
    Vec3 wall_velocity;
};

// 3D Lattice Boltzmann Grid hosting D3Q19 state and driving ISPC SPMD kernels
class LatticeGrid3D {
public:
    LatticeGrid3D(int nx, int ny, int nz, real_t relaxation_time = 0.6f);
    ~LatticeGrid3D() = default;

    LatticeGrid3D(const LatticeGrid3D&) = delete;
    LatticeGrid3D& operator=(const LatticeGrid3D&) = delete;
    LatticeGrid3D(LatticeGrid3D&&) noexcept = default;
    LatticeGrid3D& operator=(LatticeGrid3D&&) noexcept = default;

    void initialize(real_t rho_init = RHO_0, real_t ux_init = 0.0f, real_t uy_init = 0.0f, real_t uz_init = 0.0f);
    void initialize_custom(const std::function<void(int x, int y, int z, real_t& rho, real_t& ux, real_t& uy, real_t& uz)>& init_fn);

    // Single time step execution: Collision -> Streaming -> Boundaries
    void step(
        CollisionType collision = CollisionType::MRT,
        real_t gx = 0.0f,
        real_t gy = 0.0f,
        real_t gz = 0.0f,
        bool use_les = false,
        real_t smag_const = DEFAULT_SMAGORINSKY_CONSTANT
    );

    void update_macroscopic();

    // Boundary geometry setters
    void set_cell_type(int x, int y, int z, CellType type);
    CellType get_cell_type(int x, int y, int z) const;

    void set_solid_box(int x0, int y0, int z0, int x1, int y1, int z1);
    void set_solid_sphere(real_t cx, real_t cy, real_t cz, real_t radius);
    void set_solid_cylinder_z(real_t cx, real_t cy, real_t radius);
    void set_moving_wall_top(real_t u_wall_x, real_t u_wall_z = 0.0f);

    // Build curved boundary interpolation links for spherical or cylinder geometry
    void build_curved_sphere_links(real_t cx, real_t cy, real_t cz, real_t radius);

    // Accessors
    int nx() const noexcept { return nx_; }
    int ny() const noexcept { return ny_; }
    int nz() const noexcept { return nz_; }
    std::size_t total_nodes() const noexcept { return total_nodes_; }
    uint64_t current_step() const noexcept { return step_count_; }

    real_t relaxation_time() const noexcept { return tau_; }
    void set_relaxation_time(real_t tau) noexcept { tau_ = tau; }
    real_t kinematic_viscosity() const noexcept { return CS2 * (tau_ - 0.5f); }

    const real_t* rho_data() const noexcept { return rho_.data(); }
    const real_t* ux_data() const noexcept { return ux_.data(); }
    const real_t* uy_data() const noexcept { return uy_.data(); }
    const real_t* uz_data() const noexcept { return uz_.data(); }
    const real_t* speed_data() const noexcept { return speed_.data(); }
    const uint8_t* cell_type_data() const noexcept { return cell_type_.data(); }
    const real_t* f_src_data() const noexcept { return f_src_.data(); }

    real_t rho(int x, int y, int z) const { return rho_[(z * ny_ + y) * nx_ + x]; }
    real_t ux(int x, int y, int z) const { return ux_[(z * ny_ + y) * nx_ + x]; }
    real_t uy(int x, int y, int z) const { return uy_[(z * ny_ + y) * nx_ + x]; }
    real_t uz(int x, int y, int z) const { return uz_[(z * ny_ + y) * nx_ + x]; }
    real_t speed(int x, int y, int z) const { return speed_[(z * ny_ + y) * nx_ + x]; }

    real_t compute_total_mass() const;
    real_t compute_total_kinetic_energy() const;
    real_t compute_max_speed() const;

private:
    int nx_{0};
    int ny_{0};
    int nz_{0};
    std::size_t total_nodes_{0};
    uint64_t step_count_{0};
    real_t tau_{0.6f};

    bool has_moving_wall_top_{false};
    real_t u_wall_top_x_{0.0f};
    real_t u_wall_top_z_{0.0f};

    AlignedBuffer<real_t> f_src_;
    AlignedBuffer<real_t> f_dst_;

    AlignedBuffer<real_t> rho_;
    AlignedBuffer<real_t> ux_;
    AlignedBuffer<real_t> uy_;
    AlignedBuffer<real_t> uz_;
    AlignedBuffer<real_t> speed_;
    AlignedBuffer<real_t> tau_eff_;
    AlignedBuffer<uint8_t> cell_type_;

    // Curved boundary buffers for BFL kernel
    bool use_curved_bfl_{false};
    AlignedBuffer<int> bfl_fluid_indices_;
    AlignedBuffer<int> bfl_q_dirs_;
    AlignedBuffer<real_t> bfl_delta_fractions_;
    AlignedBuffer<real_t> bfl_uw_x_;
    AlignedBuffer<real_t> bfl_uw_y_;
    AlignedBuffer<real_t> bfl_uw_z_;
    int num_bfl_links_{0};
};

} // namespace wake
