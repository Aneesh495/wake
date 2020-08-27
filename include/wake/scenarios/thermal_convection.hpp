#pragma once

#include "wake/scenarios/scenario.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/thermal/thermal_solver.hpp"
#include "wake/export/terminal_renderer.hpp"
#include "wake/export/csv_telemetry.hpp"
#include <chrono>

namespace wake {

class ThermalConvectionScenario : public Scenario {
public:
    explicit ThermalConvectionScenario(
        int nx = 200,
        int ny = 100,
        real_t rayleigh = 10000.0f,
        real_t prandtl = 0.71f
    );

    ~ThermalConvectionScenario() override = default;

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
    const ThermalSolver2D& thermal() const noexcept { return thermal_solver_; }
    real_t current_nusselt() const { return thermal_solver_.compute_nusselt_number(); }

private:
    ScenarioMetadata metadata_;
    int nx_{200};
    int ny_{100};
    real_t rayleigh_{10000.0f};
    real_t prandtl_{0.71f};
    real_t g_gravity_{-0.001f};

    LatticeGrid2D grid_;
    ThermalSolver2D thermal_solver_;
    TerminalRenderer terminal_renderer_;
    TelemetryExporter telemetry_exporter_;

    double current_mlups_{0.0};
    std::chrono::high_resolution_clock::time_point last_step_time_;
};

} // namespace wake
