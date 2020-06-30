#include "wake/ibm/immersed_boundary.hpp"
#include "wake/core/constants.hpp"
#include <cmath>

namespace ispc {
extern "C" {
    void ibm_interpolate_velocity_2d(
        const float ux[], const float uy[],
        const float marker_x[], const float marker_y[],
        float marker_ux[], float marker_uy[],
        int nx, int ny, int num_markers);
    void ibm_spread_force_2d(
        const float marker_x[], const float marker_y[],
        const float marker_fx[], const float marker_fy[],
        const float marker_ds[],
        float euler_gx[], float euler_gy[],
        int nx, int ny, int num_markers);
}
}

namespace wake {

ImmersedBoundarySolver::ImmersedBoundarySolver(int nx, int ny)
    : nx_(nx), ny_(ny) {}

void ImmersedBoundarySolver::add_circle(real_t cx, real_t cy, real_t radius, int num_points) {
    if (num_points < 3) return;
    real_t d_theta = TWO_PI / static_cast<real_t>(num_points);
    real_t ds = radius * d_theta;

    for (int i = 0; i < num_points; ++i) {
        real_t theta = static_cast<real_t>(i) * d_theta;
        LagrangianMarker m;
        m.pos = Vec3(cx + radius * std::cos(theta), cy + radius * std::sin(theta), 0.0f);
        m.target_vel = Vec3(0.0f, 0.0f, 0.0f);
        m.ds = ds;
        markers_.push_back(m);
    }
    sync_buffers();
}

void ImmersedBoundarySolver::add_polygon(const std::vector<Vec2>& points) {
    if (points.size() < 2) return;
    for (std::size_t i = 0; i < points.size(); ++i) {
        std::size_t next = (i + 1) % points.size();
        Vec2 p1 = points[i];
        Vec2 p2 = points[next];
        real_t segment_len = std::sqrt((p2 - p1).norm_sq());

        LagrangianMarker m;
        m.pos = Vec3(p1.x, p1.y, 0.0f);
        m.target_vel = Vec3(0.0f, 0.0f, 0.0f);
        m.ds = segment_len;
        markers_.push_back(m);
    }
    sync_buffers();
}

void ImmersedBoundarySolver::sync_buffers() {
    std::size_t n = markers_.size();
    marker_x_.resize(n);
    marker_y_.resize(n);
    marker_ux_.resize(n);
    marker_uy_.resize(n);
    marker_fx_.resize(n);
    marker_fy_.resize(n);
    marker_ds_.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        marker_x_[i] = markers_[i].pos.x;
        marker_y_[i] = markers_[i].pos.y;
        marker_ds_[i] = markers_[i].ds;
    }
}

void ImmersedBoundarySolver::update(
    const real_t* euler_ux,
    const real_t* euler_uy,
    real_t* euler_gx,
    real_t* euler_gy,
    real_t dt
) {
    int n = static_cast<int>(markers_.size());
    if (n == 0) return;

    // 1. Interpolate Eulerian fluid velocity to Lagrangian markers via ISPC
    ispc::ibm_interpolate_velocity_2d(
        euler_ux, euler_uy,
        marker_x_.data(), marker_y_.data(),
        marker_ux_.data(), marker_uy_.data(),
        nx_, ny_, n
    );

    // 2. Compute restoring boundary force: F = (U_target - U_interp) / dt
    real_t inv_dt = 1.0f / dt;
    for (int i = 0; i < n; ++i) {
        markers_[i].fluid_vel.x = marker_ux_[i];
        markers_[i].fluid_vel.y = marker_uy_[i];

        real_t fx = (markers_[i].target_vel.x - marker_ux_[i]) * inv_dt;
        real_t fy = (markers_[i].target_vel.y - marker_uy_[i]) * inv_dt;

        marker_fx_[i] = fx;
        marker_fy_[i] = fy;

        markers_[i].force.x = fx * markers_[i].ds;
        markers_[i].force.y = fy * markers_[i].ds;
    }

    // 3. Spread Lagrangian restoring force to Eulerian lattice body force via ISPC
    ispc::ibm_spread_force_2d(
        marker_x_.data(), marker_y_.data(),
        marker_fx_.data(), marker_fy_.data(),
        marker_ds_.data(),
        euler_gx, euler_gy,
        nx_, ny_, n
    );
}

Vec3 ImmersedBoundarySolver::compute_net_force() const {
    Vec3 net_f(0.0f, 0.0f, 0.0f);
    for (const auto& m : markers_) {
        // By Newton's third law, force on body is opposite to force on fluid: F_body = -F_fluid
        net_f.x -= m.force.x;
        net_f.y -= m.force.y;
    }
    return net_f;
}

} // namespace wake
