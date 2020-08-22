#pragma once

#include "wake/scenarios/scenario.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include "wake/export/terminal_renderer.hpp"
#include "wake/export/csv_telemetry.hpp"
#include <chrono>
#include <vector>

namespace wake {

class Sphere3DWakeScenario : public Scenario {
public:
    explicit Sphere3DWakeScenario(
        int nx = 120,
        int ny = 64,
        int nz = 64,
        real_t reynolds = 300.0f,
        real_t u_inflow = 0.05f,
        real_t sphere_radius = 8.0f
    );

    ~Sphere3DWakeScenario() override = default;

    const ScenarioMetadata& metadata() const noexcept override { return metadata_; }

    void setup() override;
    void step() override;

    uint64_t current_step() const noexcept override { return grid_.current_step(); }
    real_t current_time() const noexcept override { return static_cast<real_t>(grid_.current_step()); }
    real_t current_mlups() const noexcept override { return static_cast<real_t>(current_mlups_); }

    void render_terminal(std::ostream& out) override;
    void export_frame(const std::string& output_dir, int frame_idx) override;
    void export_telemetry(const std::string& output_dir) override;

    bool is_finished() const noexcept override { return grid_.current_step() >= metadata_.default_steps; }

    const LatticeGrid3D& grid() const noexcept { return grid_; }
    real_t latest_cd() const noexcept { return latest_forces_.cd; }
    real_t latest_cl() const noexcept { return latest_forces_.cl; }

private:
    ScenarioMetadata metadata_;
    int nx_{120};
    int ny_{64};
    int nz_{64};
    real_t reynolds_{300.0f};
    real_t u_inflow_{0.05f};
    real_t sphere_radius_{8.0f};
    real_t cx_{30.0f};
    real_t cy_{32.0f};
    real_t cz_{32.0f};
    real_t tau_{0.508f};

    LatticeGrid3D grid_;
    AlignedBuffer<real_t> q_criterion_;
    AlignedBuffer<real_t> omega_x_;
    AlignedBuffer<real_t> omega_y_;
    AlignedBuffer<real_t> omega_z_;

    AerodynamicForces latest_forces_{};
    TerminalRenderer terminal_renderer_;
    TelemetryExporter telemetry_exporter_;

    double current_mlups_{0.0};
    std::chrono::high_resolution_clock::time_point last_step_time_;
};

} // namespace wake
