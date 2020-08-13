#pragma once

#include "wake/scenarios/scenario.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/geometry/naca_airfoil.hpp"
#include "wake/forces/momentum_exchange.hpp"
#include "wake/export/terminal_renderer.hpp"
#include "wake/export/csv_telemetry.hpp"
#include <chrono>
#include <vector>

namespace wake {

class NacaAirfoilScenario : public Scenario {
public:
    explicit NacaAirfoilScenario(
        const std::string& code = "0012",
        real_t angle_of_attack_deg = 10.0f,
        int nx = 480,
        int ny = 160,
        real_t reynolds = 1000.0f,
        real_t u_inflow = 0.06f,
        real_t chord = 90.0f
    );

    ~NacaAirfoilScenario() override = default;

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
    real_t latest_cm() const noexcept { return latest_forces_.cm; }

private:
    ScenarioMetadata metadata_;
    std::string airfoil_code_{"0012"};
    real_t aoa_deg_{10.0f};
    int nx_{480};
    int ny_{160};
    real_t reynolds_{1000.0f};
    real_t u_inflow_{0.06f};
    real_t chord_{90.0f};
    real_t cx_{100.0f};
    real_t cy_{80.0f};
    real_t tau_{0.5162f};

    LatticeGrid2D grid_;
    AlignedBuffer<real_t> vorticity_;
    NacaAirfoil airfoil_;

    AerodynamicForces latest_forces_{};
    TerminalRenderer terminal_renderer_;
    TelemetryExporter telemetry_exporter_;

    double current_mlups_{0.0};
    std::chrono::high_resolution_clock::time_point last_step_time_;
};

} // namespace wake
