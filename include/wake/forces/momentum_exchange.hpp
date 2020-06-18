#pragma once

#include "wake/core/types.hpp"
#include "wake/core/constants.hpp"

namespace wake {

// Aerodynamic force result record
struct AerodynamicForces {
    real_t fx{0.0f};          // Drag force (lattice units)
    real_t fy{0.0f};          // Lift force (lattice units)
    real_t fz{0.0f};          // Lateral force (3D)
    real_t torque_z{0.0f};    // Pitching moment
    real_t cd{0.0f};          // Drag coefficient
    real_t cl{0.0f};          // Lift coefficient
    real_t cm{0.0f};          // Pitching moment coefficient
};

// Calculates non-dimensional aerodynamic coefficients
inline AerodynamicForces compute_aerodynamic_coefficients(
    real_t fx, real_t fy, real_t fz, real_t torque_z,
    real_t rho_inf, real_t u_inf, real_t characteristic_area, real_t chord_length = 1.0f
) {
    AerodynamicForces forces;
    forces.fx = fx;
    forces.fy = fy;
    forces.fz = fz;
    forces.torque_z = torque_z;

    real_t dyn_pressure = 0.5f * rho_inf * u_inf * u_inf;
    if (dyn_pressure > 1e-8f && characteristic_area > 1e-8f) {
        forces.cd = fx / (dyn_pressure * characteristic_area);
        forces.cl = fy / (dyn_pressure * characteristic_area);
        forces.cm = torque_z / (dyn_pressure * characteristic_area * chord_length);
    }
    return forces;
}

} // namespace wake
