#pragma once

#include "wake/core/types.hpp"
#include "wake/core/constants.hpp"
#include "wake/core/memory.hpp"
#include "wake/grid/lattice_grid_2d.hpp"

namespace wake {

// Coupled Thermal Lattice Boltzmann solver for natural and forced convection
class ThermalSolver2D {
public:
    ThermalSolver2D(int nx, int ny, real_t tau_thermal = 0.8f, real_t t_ref = 0.5f, real_t beta = 0.05f);

    void initialize(real_t t_init = 0.5f);

    // Single thermal step: collide -> stream -> buoyancy force coupling
    void step(
        LatticeGrid2D& fluid_grid,
        real_t g_gravity = -0.001f
    );

    void set_isothermal_wall_south(real_t t_hot);
    void set_isothermal_wall_north(real_t t_cold);

    const real_t* temperature_data() const noexcept { return temperature_.data(); }
    real_t temperature(int x, int y) const { return temperature_[y * nx_ + x]; }

    real_t compute_mean_temperature() const;
    real_t compute_nusselt_number() const;

private:
    int nx_{0};
    int ny_{0};
    std::size_t total_nodes_{0};
    real_t tau_T_{0.8f};
    real_t t_ref_{0.5f};
    real_t beta_{0.05f};

    bool has_hot_south_{false};
    real_t t_hot_{1.0f};
    bool has_cold_north_{false};
    real_t t_cold_{0.0f};

    AlignedBuffer<real_t> g_src_;
    AlignedBuffer<real_t> g_dst_;
    AlignedBuffer<real_t> temperature_;
    AlignedBuffer<real_t> gy_buoyancy_;
};

} // namespace wake
