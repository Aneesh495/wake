#pragma once

#include "wake/scenarios/scenario.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/forces/momentum_exchange.hpp"
#include "wake/forces/strouhal_analyzer.hpp"
#include "wake/export/terminal_renderer.hpp"
#include "wake/export/csv_telemetry.hpp"
#include <chrono>
#include <vector>

namespace wake {

class KarmanVortexScenario : public Scenario {
public:
    explicit KarmanVortexScenario(
        int nx = 400,
        int ny = 120,
        real_t reynolds = 100.0f,
        real_t u_inflow = 0.05f,
        real_t cylinder_diameter = 20.0f
    );

    ~KarmanVortexScenario() override = default;

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

    const LatticeGrid2D& grid() const noexcept { return grid_; }
    real_t latest_cd() const noexcept { return latest_forces_.cd; }
    real_t latest_cl() const noexcept { return latest_forces_.cl; }
    real_t compute_strouhal_number() const;

private:
    ScenarioMetadata metadata_;
    int nx_{400};
    int ny_{120};
    real_t reynolds_{100.0f};
    real_t u_inflow_{0.05f};
    real_t cylinder_diameter_{20.0f};
    real_t cx_{80.0f};
    real_t cy_{60.0f};
    real_t tau_{0.53f};

    LatticeGrid2D grid_;
    AlignedBuffer<real_t> vorticity_;
    std::vector<real_t> cl_history_;

    AerodynamicForces latest_forces_{};
    mutable StrouhalAnalyzer strouhal_analyzer_;
    TerminalRenderer terminal_renderer_;
    TelemetryExporter telemetry_exporter_;

    double current_mlups_{0.0};
    std::chrono::high_resolution_clock::time_point last_step_time_;
};

} // namespace wake
