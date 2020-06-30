#pragma once

#include "wake/ibm/lagrangian_marker.hpp"
#include "wake/core/memory.hpp"
#include <vector>

namespace wake {

// Immersed Boundary Method (IBM) solver driving Lagrangian marker kinematics
class ImmersedBoundarySolver {
public:
    ImmersedBoundarySolver(int nx, int ny);

    // Add closed circle / cylinder boundary of markers
    void add_circle(real_t cx, real_t cy, real_t radius, int num_points);

    // Add arbitrary polygon or airfoil curve of markers
    void add_polygon(const std::vector<Vec2>& points);

    // Compute velocity interpolation from Eulerian lattice to markers,
    // evaluate penalty restoring force, and spread to body force field (gx, gy)
    void update(
        const real_t* euler_ux,
        const real_t* euler_uy,
        real_t* euler_gx,
        real_t* euler_gy,
        real_t dt = 1.0f
    );

    // Calculate total net aerodynamic hydrodynamic force integrated across all markers
    Vec3 compute_net_force() const;

    std::size_t num_markers() const noexcept { return markers_.size(); }
    const std::vector<LagrangianMarker>& markers() const noexcept { return markers_; }

private:
    void sync_buffers();

    int nx_{0};
    int ny_{0};
    std::vector<LagrangianMarker> markers_;

    // SIMD buffers for ISPC kernel calls
    AlignedBuffer<real_t> marker_x_;
    AlignedBuffer<real_t> marker_y_;
    AlignedBuffer<real_t> marker_ux_;
    AlignedBuffer<real_t> marker_uy_;
    AlignedBuffer<real_t> marker_fx_;
    AlignedBuffer<real_t> marker_fy_;
    AlignedBuffer<real_t> marker_ds_;
};

} // namespace wake
