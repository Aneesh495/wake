#include "wake/scenarios/thermal_convection.hpp"
#include "wake/scenarios/scenario_registry.hpp"
#include "wake/export/vtk_exporter.hpp"
#include "wake/export/image_exporter.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace wake {

ThermalConvectionScenario::ThermalConvectionScenario(
    int nx,
    int ny,
    real_t rayleigh,
    real_t prandtl
) : nx_(nx),
    ny_(ny),
    rayleigh_(rayleigh),
    prandtl_(prandtl),
    grid_(nx, ny, 0.02f),
    thermal_solver_(nx, ny, 0.8f, 0.5f, 0.05f),
    terminal_renderer_(80, 24, Colormap::Inferno) {

    metadata_.name = "rayleigh_benard_convection";
    metadata_.description = "2D natural thermal convection cells driven by Boussinesq buoyancy forces between hot bottom and cold top plates";
    metadata_.dimension = SimulationDimension::Dim2D;
    metadata_.reynolds_number = std::sqrt(rayleigh_ / prandtl_);
    metadata_.mach_number = 0.04f;
    metadata_.default_steps = 15000;
    metadata_.export_interval = 250;

    // Determine relaxation times for fluid and thermal from Pr and Ra
    real_t h = static_cast<real_t>(ny_ - 1);
    real_t delta_t = 1.0f;
    real_t beta = 0.05f;
    g_gravity_ = -0.0005f;

    // Ra = (g * beta * delta_T * H^3) / (nu * alpha)
    // Pr = nu / alpha
    real_t numerator = std::abs(g_gravity_) * beta * delta_t * (h * h * h);
    real_t nu_alpha = numerator / rayleigh_;
    real_t nu = std::sqrt(nu_alpha * prandtl_);
    real_t alpha = nu / prandtl_;

    real_t tau_fluid = 3.0f * nu + 0.5f;
    real_t tau_thermal = 3.0f * alpha + 0.5f;

    grid_.set_relaxation_time(tau_fluid);
    thermal_solver_ = ThermalSolver2D(nx, ny, tau_thermal, 0.5f, beta);
}

void ThermalConvectionScenario::setup() {
    grid_.initialize(RHO_0, 0.0f, 0.0f);

    // Initial temperature with linear conduction gradient plus small sinusoidal disturbance
    thermal_solver_.initialize(0.5f);
    thermal_solver_.set_isothermal_wall_south(1.0f); // hot bottom
    thermal_solver_.set_isothermal_wall_north(0.0f); // cold top

    last_step_time_ = std::chrono::high_resolution_clock::now();
}

void ThermalConvectionScenario::step() {
    auto t0 = std::chrono::high_resolution_clock::now();

    // 1. Step thermal solver and compute buoyancy acceleration
    thermal_solver_.step(grid_, g_gravity_);

    // 2. Step fluid solver with TRT collision coupled to vertical buoyancy force
    grid_.step(CollisionType::TRT, 0.0f, -0.0005f * 0.05f * (thermal_solver_.compute_mean_temperature() - 0.5f));

    auto t1 = std::chrono::high_resolution_clock::now();
    double step_duration_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (step_duration_ms > 1e-6) {
        current_mlups_ = (static_cast<double>(nx_ * ny_) / 1e6) / (step_duration_ms * 1e-3);
    }

    uint64_t cur_step = grid_.current_step();

    if (cur_step % 20 == 0) {
        PerformanceMetric perf;
        perf.step = cur_step;
        perf.step_time_ms = step_duration_ms;
        perf.mlups = current_mlups_;
        perf.max_speed = grid_.compute_max_speed();
        perf.total_mass = grid_.compute_total_mass();
        perf.kinetic_energy = grid_.compute_total_kinetic_energy();
        telemetry_exporter_.record_performance(perf);
    }
}

void ThermalConvectionScenario::render_terminal(std::ostream& out) {
    std::string frame = terminal_renderer_.render_2d(
        grid_,
        thermal_solver_.temperature_data(),
        0.0f,
        1.0f
    );

    std::ostringstream status;
    status << "\033[1;37m[WAKE CFD]\033[0m Scenario: Rayleigh-Benard Convection (Ra: "
           << std::fixed << std::setprecision(0) << rayleigh_ << ") | "
           << "Step: " << grid_.current_step() << " | "
           << "Throughput: " << std::setprecision(1) << current_mlups_ << " MLUPS | "
           << "Nusselt: " << std::setprecision(3) << thermal_solver_.compute_nusselt_number() << " | "
           << "Mean T: " << std::setprecision(3) << thermal_solver_.compute_mean_temperature() << "\n";

    terminal_renderer_.display(out, status.str() + frame);
}

void ThermalConvectionScenario::export_frame(const std::string& output_dir, int frame_idx) {
    std::ostringstream vti_fn;
    vti_fn << output_dir << "/thermal_" << std::setw(6) << std::setfill('0') << frame_idx << ".vti";
    VtkExporter::export_vti_2d(vti_fn.str(), grid_);

    std::ostringstream ppm_fn;
    ppm_fn << output_dir << "/thermal_temperature_" << std::setw(6) << std::setfill('0') << frame_idx << ".ppm";
    ImageExporter::export_ppm_2d(ppm_fn.str(), grid_, thermal_solver_.temperature_data(), Colormap::Inferno, 0.0f, 1.0f);
}

void ThermalConvectionScenario::export_telemetry(const std::string& output_dir) {
    telemetry_exporter_.write_performance_csv(output_dir + "/thermal_performance.csv");
    telemetry_exporter_.write_summary_json(output_dir + "/thermal_summary.json", metadata_.name, 0.0f, 0.0f);
}

REGISTER_SCENARIO("thermal", ThermalConvectionScenario)

} // namespace wake
