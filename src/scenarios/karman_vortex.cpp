#include "wake/scenarios/karman_vortex.hpp"
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

KarmanVortexScenario::KarmanVortexScenario(
    int nx,
    int ny,
    real_t reynolds,
    real_t u_inflow,
    real_t cylinder_diameter
) : nx_(nx),
    ny_(ny),
    reynolds_(reynolds),
    u_inflow_(u_inflow),
    cylinder_diameter_(cylinder_diameter),
    grid_(nx, ny, 0.02f),
    vorticity_(nx * ny),
    terminal_renderer_(80, 24, Colormap::CoolWarm) {

    metadata_.name = "karman_vortex_street";
    metadata_.description = "2D unsteady flow past a circular cylinder with laminar periodic Karman vortex shedding";
    metadata_.dimension = SimulationDimension::Dim2D;
    metadata_.reynolds_number = reynolds_;
    metadata_.mach_number = u_inflow_ / CS;
    metadata_.default_steps = 15000;
    metadata_.export_interval = 200;

    cx_ = static_cast<real_t>(nx_) * 0.20f;
    cy_ = static_cast<real_t>(ny_) * 0.50f;

    // Kinematic viscosity: nu = u * D / Re
    real_t nu = (u_inflow_ * cylinder_diameter_) / reynolds_;
    // Relaxation time: tau = 3 * nu + 0.5
    tau_ = 3.0f * nu + 0.5f;
    grid_.set_relaxation_time(tau_);
}

void KarmanVortexScenario::setup() {
    grid_.set_relaxation_time(tau_);

    // Initialize with freestream flow plus small asymmetry to trigger shedding
    grid_.initialize_custom([this](int x, int y, real_t& rho, real_t& ux, real_t& uy) {
        rho = RHO_0;
        ux = u_inflow_;
        real_t dx = static_cast<real_t>(x) - cx_;
        real_t dy = static_cast<real_t>(y) - cy_;
        real_t r2 = dx * dx + dy * dy;
        real_t perturb = 0.005f * u_inflow_ * std::sin(PI * dy / static_cast<real_t>(ny_));
        uy = (r2 < cylinder_diameter_ * cylinder_diameter_ * 4.0f) ? perturb : 0.0f;
    });

    // Inflow on West (x = 0), Outflow on East (x = nx - 1)
    grid_.set_inflow_west(u_inflow_, 0.0f);
    grid_.set_outflow_east(RHO_0);

    // Rasterize circular cylinder obstacle
    real_t radius = cylinder_diameter_ * 0.5f;
    grid_.set_solid_circle(cx_, cy_, radius);

    cl_history_.clear();
    last_step_time_ = std::chrono::high_resolution_clock::now();
}

void KarmanVortexScenario::step() {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Execute one LBM time step using Two-Relaxation-Time (TRT) collision
    grid_.step(CollisionType::TRT);

    auto t1 = std::chrono::high_resolution_clock::now();
    double step_duration_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (step_duration_ms > 1e-6) {
        current_mlups_ = (static_cast<double>(nx_ * ny_) / 1e6) / (step_duration_ms * 1e-3);
    }

    uint64_t cur_step = grid_.current_step();

    // Compute hydrodynamic forces via Momentum Exchange Method
    float fx = 0.0f, fy = 0.0f, tz = 0.0f;
    ispc::d2q9_compute_momentum_exchange_forces(
        grid_.f_src_data(),
        grid_.cell_type_data(),
        nx_, ny_,
        cx_, cy_,
        fx, fy, tz
    );

    latest_forces_ = compute_aerodynamic_coefficients(
        fx, fy, 0.0f, tz,
        RHO_0, u_inflow_,
        cylinder_diameter_, cylinder_diameter_
    );

    // Record lift coefficient sample for Strouhal analysis
    cl_history_.push_back(latest_forces_.cl);

    // Record telemetry periodically
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

real_t KarmanVortexScenario::compute_strouhal_number() const {
    if (cl_history_.size() < 64) {
        return 0.0f;
    }
    real_t f = strouhal_analyzer_.analyze_frequency(cl_history_, 1.0f);
    return strouhal_analyzer_.compute_strouhal_number(f, cylinder_diameter_, u_inflow_);
}

void KarmanVortexScenario::render_terminal(std::ostream& out) {
    FieldAnalytics::compute_vorticity(grid_, vorticity_);
    std::string frame = terminal_renderer_.render_vorticity_2d(grid_, vorticity_.data(), 0.03f);

    std::ostringstream status;
    status << "\033[1;37m[WAKE CFD]\033[0m Scenario: Karman Vortex Street | "
           << "Step: " << grid_.current_step() << " | "
           << "Throughput: " << std::fixed << std::setprecision(1) << current_mlups_ << " MLUPS | "
           << "Cd: " << std::setprecision(3) << latest_forces_.cd << " | "
           << "Cl: " << std::setprecision(3) << latest_forces_.cl << "\n";

    terminal_renderer_.display(out, status.str() + frame);
}

void KarmanVortexScenario::export_frame(const std::string& output_dir, int frame_idx) {
    FieldAnalytics::compute_vorticity(grid_, vorticity_);

    std::ostringstream vti_fn;
    vti_fn << output_dir << "/karman_" << std::setw(6) << std::setfill('0') << frame_idx << ".vti";
    VtkExporter::export_vti_2d(vti_fn.str(), grid_, vorticity_.data());

    std::ostringstream ppm_fn;
    ppm_fn << output_dir << "/karman_vorticity_" << std::setw(6) << std::setfill('0') << frame_idx << ".ppm";
    ImageExporter::export_ppm_vorticity_2d(ppm_fn.str(), grid_, vorticity_.data(), 0.03f);
}

void KarmanVortexScenario::export_telemetry(const std::string& output_dir) {
    telemetry_exporter_.write_aero_csv(output_dir + "/karman_aerodynamics.csv");
    telemetry_exporter_.write_performance_csv(output_dir + "/karman_performance.csv");

    real_t f = 0.0f;
    real_t st = 0.0f;
    if (cl_history_.size() >= 64) {
        f = strouhal_analyzer_.analyze_frequency(cl_history_, 1.0f);
        st = strouhal_analyzer_.compute_strouhal_number(f, cylinder_diameter_, u_inflow_);
    }
    telemetry_exporter_.write_summary_json(output_dir + "/karman_summary.json", metadata_.name, st, f);
}

REGISTER_SCENARIO("karman", KarmanVortexScenario)

} // namespace wake
