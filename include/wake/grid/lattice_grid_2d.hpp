#pragma once

#include "wake/core/types.hpp"
#include "wake/core/constants.hpp"
#include "wake/core/memory.hpp"
#include "wake/core/lattice_models.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace wake {

// 2D Lattice Boltzmann Grid hosting D2Q9 state and driving ISPC SPMD kernels
class LatticeGrid2D {
public:
    LatticeGrid2D(int nx, int ny, real_t relaxation_time = 0.6f);
    ~LatticeGrid2D() = default;

    // Non-copyable, movable
    LatticeGrid2D(const LatticeGrid2D&) = delete;
    LatticeGrid2D& operator=(const LatticeGrid2D&) = delete;
    LatticeGrid2D(LatticeGrid2D&&) noexcept = default;
    LatticeGrid2D& operator=(LatticeGrid2D&&) noexcept = default;

    // Initialization with uniform or prescribed macroscopic field
    void initialize(real_t rho_init = RHO_0, real_t ux_init = 0.0f, real_t uy_init = 0.0f);

    // Advanced initialization with custom initializer functor
    void initialize_custom(const std::function<void(int x, int y, real_t& rho, real_t& ux, real_t& uy)>& init_fn);

    // Single time step execution: Collision -> Streaming -> Boundary Conditions
    void step(
        CollisionType collision = CollisionType::MRT,
        real_t gx = 0.0f,
        real_t gy = 0.0f,
        bool use_les = false,
        real_t smag_const = DEFAULT_SMAGORINSKY_CONSTANT
    );

    // Update macroscopic density and velocity from current distributions
    void update_macroscopic();

    // Boundary configuration methods
    void set_cell_type(int x, int y, CellType type);
    CellType get_cell_type(int x, int y) const;

    void set_solid_box(int x0, int y0, int x1, int y1);
    void set_solid_circle(real_t cx, real_t cy, real_t radius);
    void set_inflow_west(real_t u_in_x, real_t u_in_y = 0.0f);
    void set_outflow_east(real_t rho_out = RHO_0);
    void set_moving_wall_north(real_t u_wall_x);

    // Accessors for physical fields
    int nx() const noexcept { return nx_; }
    int ny() const noexcept { return ny_; }
    std::size_t total_nodes() const noexcept { return total_nodes_; }
    uint64_t current_step() const noexcept { return step_count_; }

    real_t relaxation_time() const noexcept { return tau_; }
    void set_relaxation_time(real_t tau) noexcept { tau_ = tau; }

    real_t kinematic_viscosity() const noexcept { return CS2 * (tau_ - 0.5f); }

    // Direct buffer pointers for high-throughput reads and visualization
    const real_t* rho_data() const noexcept { return rho_.data(); }
    const real_t* ux_data() const noexcept { return ux_.data(); }
    const real_t* uy_data() const noexcept { return uy_.data(); }
    const real_t* speed_data() const noexcept { return speed_.data(); }
    const uint8_t* cell_type_data() const noexcept { return cell_type_.data(); }
    const real_t* f_src_data() const noexcept { return f_src_.data(); }
    const real_t* f_dst_data() const noexcept { return f_dst_.data(); }

    real_t rho(int x, int y) const { return rho_[y * nx_ + x]; }
    real_t ux(int x, int y) const { return ux_[y * nx_ + x]; }
    real_t uy(int x, int y) const { return uy_[y * nx_ + x]; }
    real_t speed(int x, int y) const { return speed_[y * nx_ + x]; }

    // Physical conservation diagnostics
    real_t compute_total_mass() const;
    real_t compute_total_kinetic_energy() const;
    real_t compute_max_speed() const;

private:
    int nx_{0};
    int ny_{0};
    std::size_t total_nodes_{0};
    uint64_t step_count_{0};
    real_t tau_{0.6f};

    // Boundary configuration flags
    bool has_inflow_west_{false};
    real_t u_inflow_x_{0.0f};
    real_t u_inflow_y_{0.0f};

    bool has_outflow_east_{false};
    real_t rho_outflow_{RHO_0};

    bool has_moving_wall_north_{false};
    real_t u_wall_north_x_{0.0f};

    // Primary distribution functions (SoA ping-pong layout)
    AlignedBuffer<real_t> f_src_;
    AlignedBuffer<real_t> f_dst_;
    AlignedBuffer<real_t> f_prev_outflow_;

    // Macroscopic physical fields
    AlignedBuffer<real_t> rho_;
    AlignedBuffer<real_t> ux_;
    AlignedBuffer<real_t> uy_;
    AlignedBuffer<real_t> speed_;
    AlignedBuffer<real_t> tau_eff_;
    AlignedBuffer<uint8_t> cell_type_;
};

} // namespace wake
