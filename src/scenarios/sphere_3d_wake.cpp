#include "wake/scenarios/sphere_3d_wake.hpp"
#include "wake/scenarios/scenario_registry.hpp"
#include "wake/analytics/field_analytics.hpp"
#include "wake/export/vtk_exporter.hpp"
#include "wake/export/image_exporter.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace ispc {
extern "C" {
    void d3q19_compute_momentum_exchange_forces(
        const float f_post_collision[],
        const uint8_t cell_type[],
        int nx, int ny, int nz,
        float center_x, float center_y, float center_z,
        float max_radius,
        float &fx_out, float &fy_out, float &fz_out
    );
}
}

namespace wake {

Sphere3DWakeScenario::Sphere3DWakeScenario(
    int nx,
    int ny,
    int nz,
    real_t reynolds,
    real_t u_inflow,
    real_t sphere_radius
) : nx_(nx),
    ny_(ny),
    nz_(nz),
    reynolds_(reynolds),
    u_inflow_(u_inflow),
    sphere_radius_(sphere_radius),
    grid_(nx, ny, nz, 0.02f),
    q_criterion_(nx * ny * nz),
    omega_x_(nx * ny * nz),
    omega_y_(nx * ny * nz),
    omega_z_(nx * ny * nz),
    terminal_renderer_(80, 24, Colormap::Inferno) {

    metadata_.name = "sphere_3d_turbulent_wake";
    metadata_.description = "3D flow past a sphere featuring vortex hairpin detachment and Smagorinsky subgrid turbulence";
    metadata_.dimension = SimulationDimension::Dim3D;
    metadata_.reynolds_number = reynolds_;
    metadata_.mach_number = u_inflow_ / CS;
    metadata_.default_steps = 10000;
    metadata_.export_interval = 200;

    cx_ = static_cast<real_t>(nx_) * 0.25f;
    cy_ = static_cast<real_t>(ny_) * 0.50f;
    cz_ = static_cast<real_t>(nz_) * 0.50f;

    // Characteristic diameter D = 2 * R
    real_t d = 2.0f * sphere_radius_;
    real_t nu = (u_inflow_ * d) / reynolds_;
    tau_ = 3.0f * nu + 0.5f;
    grid_.set_relaxation_time(tau_);
}

void Sphere3DWakeScenario::setup() {
    grid_.set_relaxation_time(tau_);
    grid_.set_periodic(true); // periodic spanwise in y and z

    // Initialize with freestream flow
    grid_.initialize(RHO_0, u_inflow_, 0.0f, 0.0f);

    // Voxelize solid sphere obstacle
    grid_.set_solid_sphere(cx_, cy_, cz_, sphere_radius_);

    // Build curved boundary links for BFL interpolation
    grid_.build_curved_sphere_links(cx_, cy_, cz_, sphere_radius_);

    last_step_time_ = std::chrono::high_resolution_clock::now();
}

void Sphere3DWakeScenario::step() {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Step with Smagorinsky LES turbulence modeling enabled at higher Reynolds numbers
    bool use_les = (reynolds_ >= 500.0f);
    grid_.step(CollisionType::MRT, 0.0f, 0.0f, 0.0f, use_les, 0.12f);

    auto t1 = std::chrono::high_resolution_clock::now();
    double step_duration_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (step_duration_ms > 1e-6) {
        current_mlups_ = (static_cast<double>(nx_ * ny_ * nz_) / 1e6) / (step_duration_ms * 1e-3);
    }

    uint64_t cur_step = grid_.current_step();

    // Compute 3D aerodynamic forces
    float fx = 0.0f, fy = 0.0f, fz = 0.0f;
    ispc::d3q19_compute_momentum_exchange_forces(
        grid_.f_src_data(),
        grid_.cell_type_data(),
        nx_, ny_, nz_,
        cx_, cy_, cz_,
        sphere_radius_ * 2.0f,
        fx, fy, fz
    );

    real_t d = 2.0f * sphere_radius_;
    real_t frontal_area = PI * sphere_radius_ * sphere_radius_;

    latest_forces_ = compute_aerodynamic_coefficients(
        fx, fy, fz, 0.0f,
        RHO_0, u_inflow_,
        frontal_area, d
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

void Sphere3DWakeScenario::render_terminal(std::ostream& out) {
    // Render mid-plane slice (z = nz / 2)
    std::string frame = terminal_renderer_.render_slice_3d(grid_, nz_ / 2, 0.0f, u_inflow_ * 1.5f);

    std::ostringstream status;
    status << "\033[1;37m[WAKE CFD]\033[0m Scenario: 3D Sphere Wake (Re: "
           << std::fixed << std::setprecision(0) << reynolds_ << ") | "
           << "Step: " << grid_.current_step() << " | "
           << "Throughput: " << std::setprecision(1) << current_mlups_ << " MLUPS | "
           << "Cd: " << std::setprecision(3) << latest_forces_.cd << " | "
           << "Cl: " << std::setprecision(3) << latest_forces_.cl << "\n";

    terminal_renderer_.display(out, status.str() + frame);
}

void Sphere3DWakeScenario::export_frame(const std::string& output_dir, int frame_idx) {
    // Compute Q-criterion and vorticity vector fields in parallel ISPC
    FieldAnalytics::compute_3d_features(grid_, omega_x_, omega_y_, omega_z_, q_criterion_);

    std::ostringstream vti_fn;
    vti_fn << output_dir << "/sphere_" << std::setw(6) << std::setfill('0') << frame_idx << ".vti";
    VtkExporter::export_vti_3d(vti_fn.str(), grid_, q_criterion_.data());

    std::ostringstream ppm_fn;
    ppm_fn << output_dir << "/sphere_midplane_" << std::setw(6) << std::setfill('0') << frame_idx << ".ppm";
    ImageExporter::export_ppm_slice_xy(ppm_fn.str(), grid_, nz_ / 2, Colormap::Inferno, 0.0f, u_inflow_ * 1.5f);
}

void Sphere3DWakeScenario::export_telemetry(const std::string& output_dir) {
    telemetry_exporter_.write_aero_csv(output_dir + "/sphere_aerodynamics.csv");
    telemetry_exporter_.write_performance_csv(output_dir + "/sphere_performance.csv");
    telemetry_exporter_.write_summary_json(output_dir + "/sphere_summary.json", metadata_.name, 0.0f, 0.0f);
}

REGISTER_SCENARIO("sphere", Sphere3DWakeScenario)

} // namespace wake
