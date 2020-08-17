#pragma once

#include "wake/scenarios/scenario.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include "wake/export/terminal_renderer.hpp"
#include "wake/export/csv_telemetry.hpp"
#include <chrono>
#include <vector>

namespace wake {

class LidDrivenCavityScenario : public Scenario {
public:
    explicit LidDrivenCavityScenario(
        int n = 48,
        real_t reynolds = 400.0f,
        real_t u_lid = 0.08f
    );

    ~LidDrivenCavityScenario() override = default;

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

    // Centerline velocity profiles for benchmark verification
    std::vector<real_t> extract_vertical_centerline_ux() const;
    std::vector<real_t> extract_horizontal_centerline_uy() const;

private:
    ScenarioMetadata metadata_;
    int n_{48};
    real_t reynolds_{400.0f};
    real_t u_lid_{0.08f};
    real_t tau_{0.5288f};

    LatticeGrid3D grid_;
    TerminalRenderer terminal_renderer_;
    TelemetryExporter telemetry_exporter_;

    double current_mlups_{0.0};
    std::chrono::high_resolution_clock::time_point last_step_time_;
};

} // namespace wake
