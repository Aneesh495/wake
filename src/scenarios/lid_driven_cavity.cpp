#include "wake/scenarios/lid_driven_cavity.hpp"
#include "wake/scenarios/scenario_registry.hpp"
#include "wake/export/vtk_exporter.hpp"
#include "wake/export/image_exporter.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace wake {

LidDrivenCavityScenario::LidDrivenCavityScenario(
    int n,
    real_t reynolds,
    real_t u_lid
) : n_(n),
    reynolds_(reynolds),
    u_lid_(u_lid),
    grid_(n, n, n, 0.02f),
    terminal_renderer_(80, 24, Colormap::Viridis) {

    metadata_.name = "lid_driven_cavity_3d";
    metadata_.description = "3D cubic lid-driven cavity benchmark problem comparing centerline profiles against Ghia et al.";
    metadata_.dimension = SimulationDimension::Dim3D;
    metadata_.reynolds_number = reynolds_;
    metadata_.mach_number = u_lid_ / CS;
    metadata_.default_steps = 10000;
    metadata_.export_interval = 200;

    // Kinematic viscosity: nu = u * L / Re
    real_t nu = (u_lid_ * static_cast<real_t>(n_)) / reynolds_;
    tau_ = 3.0f * nu + 0.5f;
    grid_.set_relaxation_time(tau_);
}

void LidDrivenCavityScenario::setup() {
    grid_.set_relaxation_time(tau_);
    grid_.set_periodic(false);

    // Initialize resting fluid
    grid_.initialize(RHO_0, 0.0f, 0.0f, 0.0f);

    // Enforce no-slip boundary walls on bottom, left, right, front, back
    for (int z = 0; z < n_; ++z) {
        for (int x = 0; x < n_; ++x) {
            grid_.set_cell_type(x, 0, z, CellType::Solid); // bottom
        }
    }
    for (int z = 0; z < n_; ++z) {
        for (int y = 0; y < n_ - 1; ++y) {
            grid_.set_cell_type(0, y, z, CellType::Solid);       // west
            grid_.set_cell_type(n_ - 1, y, z, CellType::Solid); // east
        }
    }
    for (int y = 0; y < n_ - 1; ++y) {
        for (int x = 0; x < n_; ++x) {
            grid_.set_cell_type(x, y, 0, CellType::Solid);       // front
            grid_.set_cell_type(x, y, n_ - 1, CellType::Solid); // back
        }
    }

    // Top moving lid
    grid_.set_moving_wall_top(u_lid_, 0.0f);

    last_step_time_ = std::chrono::high_resolution_clock::now();
}

void LidDrivenCavityScenario::step() {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Step with 3D MRT collision operator
    grid_.step(CollisionType::MRT);

    auto t1 = std::chrono::high_resolution_clock::now();
    double step_duration_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (step_duration_ms > 1e-6) {
        current_mlups_ = (static_cast<double>(n_ * n_ * n_) / 1e6) / (step_duration_ms * 1e-3);
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

std::vector<real_t> LidDrivenCavityScenario::extract_vertical_centerline_ux() const {
    int mid_x = n_ / 2;
    int mid_z = n_ / 2;
    std::vector<real_t> ux_profile(n_);
    for (int y = 0; y < n_; ++y) {
        ux_profile[y] = grid_.ux(mid_x, y, mid_z);
    }
    return ux_profile;
}

std::vector<real_t> LidDrivenCavityScenario::extract_horizontal_centerline_uy() const {
    int mid_y = n_ / 2;
    int mid_z = n_ / 2;
    std::vector<real_t> uy_profile(n_);
    for (int x = 0; x < n_; ++x) {
        uy_profile[x] = grid_.uy(x, mid_y, mid_z);
    }
    return uy_profile;
}

void LidDrivenCavityScenario::render_terminal(std::ostream& out) {
    std::string frame = terminal_renderer_.render_slice_3d(grid_, n_ / 2, 0.0f, u_lid_);

    std::ostringstream status;
    status << "\033[1;37m[WAKE CFD]\033[0m Scenario: 3D Lid-Driven Cavity (Re: "
           << std::fixed << std::setprecision(0) << reynolds_ << ") | "
           << "Step: " << grid_.current_step() << " | "
           << "Throughput: " << std::setprecision(1) << current_mlups_ << " MLUPS | "
           << "Max Speed: " << std::setprecision(4) << grid_.compute_max_speed() << "\n";

    terminal_renderer_.display(out, status.str() + frame);
}

void LidDrivenCavityScenario::export_frame(const std::string& output_dir, int frame_idx) {
    std::ostringstream vti_fn;
    vti_fn << output_dir << "/cavity_" << std::setw(6) << std::setfill('0') << frame_idx << ".vti";
    VtkExporter::export_vti_3d(vti_fn.str(), grid_);

    std::ostringstream ppm_fn;
    ppm_fn << output_dir << "/cavity_midplane_" << std::setw(6) << std::setfill('0') << frame_idx << ".ppm";
    ImageExporter::export_ppm_slice_xy(ppm_fn.str(), grid_, n_ / 2, Colormap::Viridis, 0.0f, u_lid_);
}

void LidDrivenCavityScenario::export_telemetry(const std::string& output_dir) {
    telemetry_exporter_.write_performance_csv(output_dir + "/cavity_performance.csv");
    telemetry_exporter_.write_summary_json(output_dir + "/cavity_summary.json", metadata_.name, 0.0f, 0.0f);
}

REGISTER_SCENARIO("cavity", LidDrivenCavityScenario)

} // namespace wake
