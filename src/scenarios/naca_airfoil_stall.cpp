#include "wake/scenarios/naca_airfoil_stall.hpp"
#include "wake/scenarios/scenario_registry.hpp"
#include "wake/analytics/field_analytics.hpp"
#include "wake/export/vtk_exporter.hpp"
#include "wake/export/image_exporter.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace ispc {
extern "C" {
    void d2q9_compute_momentum_exchange_forces(
        const float f_post_collision[],
        const uint8_t cell_type[],
        int nx, int ny,
        float center_x, float center_y,
        float &fx_out, float &fy_out, float &torque_z_out
    );
}
}

namespace wake {

NacaAirfoilScenario::NacaAirfoilScenario(
    const std::string& code,
    real_t angle_of_attack_deg,
    int nx,
    int ny,
    real_t reynolds,
    real_t u_inflow,
    real_t chord
) : airfoil_code_(code),
    aoa_deg_(angle_of_attack_deg),
    nx_(nx),
    ny_(ny),
    reynolds_(reynolds),
    u_inflow_(u_inflow),
    chord_(chord),
    grid_(nx, ny, 0.02f),
    vorticity_(nx * ny),
    airfoil_(code, chord, angle_of_attack_deg),
    terminal_renderer_(80, 24, Colormap::CoolWarm) {

    metadata_.name = "naca_airfoil_stall";
    metadata_.description = "2D aerodynamic flow and dynamic stall analysis past a NACA 4-digit airfoil";
    metadata_.dimension = SimulationDimension::Dim2D;
    metadata_.reynolds_number = reynolds_;
    metadata_.mach_number = u_inflow_ / CS;
    metadata_.default_steps = 20000;
    metadata_.export_interval = 250;

    cx_ = static_cast<real_t>(nx_) * 0.22f;
    cy_ = static_cast<real_t>(ny_) * 0.50f;

    // Kinematic viscosity: nu = u * chord / Re
    real_t nu = (u_inflow_ * chord_) / reynolds_;
    tau_ = 3.0f * nu + 0.5f;
    grid_.set_relaxation_time(tau_);
}

void NacaAirfoilScenario::setup() {
    grid_.set_relaxation_time(tau_);

    // Initialize with freestream flow
    grid_.initialize(RHO_0, u_inflow_, 0.0f);

    // Boundaries: West inflow, East outflow
    grid_.set_inflow_west(u_inflow_, 0.0f);
    grid_.set_outflow_east(RHO_0);

    // Rasterize airfoil geometry with angle of attack
    airfoil_.rasterize(grid_, cx_, cy_);

    last_step_time_ = std::chrono::high_resolution_clock::now();
}

void NacaAirfoilScenario::step() {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Step with TRT collision operator
    grid_.step(CollisionType::TRT);

    auto t1 = std::chrono::high_resolution_clock::now();
    double step_duration_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (step_duration_ms > 1e-6) {
        current_mlups_ = (static_cast<double>(nx_ * ny_) / 1e6) / (step_duration_ms * 1e-3);
    }

    uint64_t cur_step = grid_.current_step();

    // Quarter-chord aerodynamic center reference point
    real_t qc_x = cx_ + chord_ * 0.25f;
    real_t qc_y = cy_;

    float fx = 0.0f, fy = 0.0f, tz = 0.0f;
    ispc::d2q9_compute_momentum_exchange_forces(
        grid_.f_src_data(),
        grid_.cell_type_data(),
        nx_, ny_,
        qc_x, qc_y,
        fx, fy, tz
    );

    latest_forces_ = compute_aerodynamic_coefficients(
        fx, fy, 0.0f, tz,
        RHO_0, u_inflow_,
        chord_, chord_
    );

    if (cur_step % 10 == 0) {
        telemetry_exporter_.record_aero({
            cur_step,
            static_cast<real_t>(cur_step),
            latest_forces_
        });

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

void NacaAirfoilScenario::render_terminal(std::ostream& out) {
    FieldAnalytics::compute_vorticity(grid_, vorticity_);
    std::string frame = terminal_renderer_.render_vorticity_2d(grid_, vorticity_.data(), 0.04f);

    std::ostringstream status;
    status << "\033[1;37m[WAKE CFD]\033[0m Scenario: NACA " << airfoil_code_ << " (AoA: "
           << std::fixed << std::setprecision(1) << aoa_deg_ << " deg) | "
           << "Step: " << grid_.current_step() << " | "
           << "Throughput: " << current_mlups_ << " MLUPS | "
           << "Cd: " << std::setprecision(3) << latest_forces_.cd << " | "
           << "Cl: " << std::setprecision(3) << latest_forces_.cl << " | "
           << "Cm: " << std::setprecision(3) << latest_forces_.cm << "\n";

    terminal_renderer_.display(out, status.str() + frame);
}

void NacaAirfoilScenario::export_frame(const std::string& output_dir, int frame_idx) {
    FieldAnalytics::compute_vorticity(grid_, vorticity_);

    std::ostringstream vti_fn;
    vti_fn << output_dir << "/airfoil_" << std::setw(6) << std::setfill('0') << frame_idx << ".vti";
    VtkExporter::export_vti_2d(vti_fn.str(), grid_, vorticity_.data());

    std::ostringstream ppm_fn;
    ppm_fn << output_dir << "/airfoil_vorticity_" << std::setw(6) << std::setfill('0') << frame_idx << ".ppm";
    ImageExporter::export_ppm_vorticity_2d(ppm_fn.str(), grid_, vorticity_.data(), 0.04f);
}

void NacaAirfoilScenario::export_telemetry(const std::string& output_dir) {
    telemetry_exporter_.write_aero_csv(output_dir + "/airfoil_aerodynamics.csv");
    telemetry_exporter_.write_performance_csv(output_dir + "/airfoil_performance.csv");
    telemetry_exporter_.write_summary_json(output_dir + "/airfoil_summary.json", metadata_.name, 0.0f, 0.0f);
}

REGISTER_SCENARIO("airfoil", NacaAirfoilScenario)

} // namespace wake
