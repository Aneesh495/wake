#pragma once

#include "wake/core/types.hpp"

namespace wake {

// A single Lagrangian boundary marker point in 2D or 3D space
struct LagrangianMarker {
    Vec3 pos{0.0f, 0.0f, 0.0f};          // Current position
    Vec3 target_vel{0.0f, 0.0f, 0.0f};    // Prescribed target boundary velocity
    Vec3 fluid_vel{0.0f, 0.0f, 0.0f};     // Interpolated local fluid velocity
    Vec3 force{0.0f, 0.0f, 0.0f};         // Computed restoring force
    real_t ds{1.0f};                      // Boundary element arc length / area weighting
};

} // namespace wake
