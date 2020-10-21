#include "wake/core/constants.hpp"
#include "wake/scenarios/scenario_registry.hpp"
#include "wake/runner/simulation_runner.hpp"
#include "wake/geometry/obj_loader.hpp"
#include "wake/export/vtk_exporter.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <csignal>

namespace {

void print_banner() {
    std::cout << "\033[1;36m"
              << " __          __     _  ________ \n"
              << " \\ \\        / /\\   | |/ /  ____|\n"
              << "  \\ \\  /\\  / /  \\  | ' /| |__   \n"
              << "   \\ \\/  \\/ / /\\ \\ |  < |  __|  \n"
              << "    \\  /\\  / ____ \\| . \\| |____ \n"
              << "     \\/  \\/_/    \\_\\_|\\_\\______|\n"
              << "\033[0m"
              << "\033[1;37m High-Performance ISPC SPMD Lattice Boltzmann CFD Engine\033[0m\n"
              << " Version: 1.0.0 | Architecture: "
#if defined(__aarch64__) || defined(_M_ARM64)
              << "arm64 (neon-i32x4)"
#elif defined(__x86_64__) || defined(_M_X64)
              << "x86_64 (avx2-i32x8 / avx512skx-x16)"
#else
              << "generic"
#endif
              << " | Alignment: " << wake::SIMD_ALIGNMENT << " bytes\n\n";
}

void print_usage(const char* prog) {
    print_banner();
    std::cout << "Usage: " << prog << " <command> [options]\n\n"
              << "Commands:\n"
              << "  run <scenario>       Execute a CFD simulation scenario\n"
              << "  list                 List all registered simulation scenarios\n"
              << "  mesh <obj_file>      Voxelize a 3D Wavefront OBJ mesh into VTK ImageData (.vti)\n"
              << "  version              Display build and hardware SIMD vectorization details\n"
              << "  help                 Display this manual\n\n"
              << "Run Options:\n"
              << "  --steps <N>          Maximum simulation timesteps (default: scenario default)\n"
              << "  --export <dir>       Directory for VTK (.vti), PPM, and CSV output (default: ./output)\n"
              << "  --interval <N>       Data export frame interval (default: 200)\n"
              << "  --render-rate <N>    Terminal live rendering cadence in steps (default: 50)\n"
              << "  --headless           Disable ANSI TrueColor terminal renderer (max throughput)\n"
              << "  --fps <N>            Throttle terminal animation to target FPS (default: unthrottled)\n\n"
              << "Examples:\n"
              << "  " << prog << " run karman\n"
              << "  " << prog << " run airfoil --steps 5000 --export ./results\n"
              << "  " << prog << " run cavity --headless\n"
              << "  " << prog << " run sphere --steps 1000\n"
              << "  " << prog << " run thermal\n"
              << "  " << prog << " list\n"
              << "  " << prog << " mesh model.obj --grid 64\n\n";
}

std::shared_ptr<wake::SimulationRunner> g_active_runner = nullptr;

void signal_handler(int sig) {
    if (sig == SIGINT || sig == SIGTERM) {
        std::cout << "\n[WAKE] Received termination signal. Gracefully finalizing simulation...\n";
        if (g_active_runner) {
            g_active_runner->stop();
        }
    }
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string cmd = argv[1];

    if (cmd == "-h" || cmd == "--help" || cmd == "help") {
        print_usage(argv[0]);
        return 0;
    }

    if (cmd == "-v" || cmd == "--version" || cmd == "version") {
        print_banner();
        return 0;
    }

    if (cmd == "list") {
        print_banner();
        std::cout << "Registered Simulation Scenarios:\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        auto list = wake::ScenarioRegistry::instance().list_scenarios();
        for (const auto& name : list) {
            auto scn = wake::ScenarioRegistry::instance().create(name);
            const auto& meta = scn->metadata();
            std::cout << "  * \033[1;32m" << std::setw(12) << std::left << name << "\033[0m: "
                      << meta.description << "\n"
                      << "    [Dimension: " << (meta.dimension == wake::SimulationDimension::Dim2D ? "2D" : "3D")
                      << ", Re: " << meta.reynolds_number
                      << ", Default Steps: " << meta.default_steps << "]\n\n";
        }
        return 0;
    }

    if (cmd == "run") {
        if (argc < 3) {
            std::cerr << "Error: 'run' requires a scenario name. Use '" << argv[0] << " list' to view scenarios.\n";
            return 1;
        }

        std::string scenario_name = "";
        std::string export_dir_override = "";
        uint64_t steps_override = 0;
        uint64_t interval_override = 0;
        uint64_t render_interval_override = 50;
        bool live_render = true;
        double target_fps = 0.0;

        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            if ((arg == "--scenario" || arg == "-s") && i + 1 < argc) {
                scenario_name = argv[++i];
            } else if ((arg == "--steps" || arg == "-n") && i + 1 < argc) {
                steps_override = std::stoull(argv[++i]);
            } else if ((arg == "--export" || arg == "--output" || arg == "-o") && i + 1 < argc) {
                export_dir_override = argv[++i];
            } else if ((arg == "--interval" || arg == "-i") && i + 1 < argc) {
                interval_override = std::stoull(argv[++i]);
            } else if (arg == "--render-rate" && i + 1 < argc) {
                render_interval_override = std::stoull(argv[++i]);
            } else if (arg == "--headless") {
                live_render = false;
            } else if (arg == "--fps" && i + 1 < argc) {
                target_fps = std::stod(argv[++i]);
            } else if (scenario_name.empty() && !arg.starts_with("-")) {
                scenario_name = arg;
            }
        }

        if (scenario_name.empty()) {
            std::cerr << "Error: No scenario specified. Use '" << argv[0] << " list' to view scenarios.\n";
            return 1;
        }

        if (!wake::ScenarioRegistry::instance().has_scenario(scenario_name)) {
            std::cerr << "Error: Unknown scenario '" << scenario_name << "'. Available scenarios:\n";
            for (const auto& name : wake::ScenarioRegistry::instance().list_scenarios()) {
                std::cerr << "  - " << name << "\n";
            }
            return 1;
        }

        auto scenario = wake::ScenarioRegistry::instance().create(scenario_name);
        const auto& meta = scenario->metadata();

        wake::RunnerConfig config;
        config.max_steps = steps_override > 0 ? steps_override : meta.default_steps;
        config.export_interval = interval_override > 0 ? interval_override : meta.export_interval;
        config.render_interval = render_interval_override;
        config.live_terminal_render = live_render;
        config.export_dir = !export_dir_override.empty() ? export_dir_override : ("./output_" + scenario_name);
        config.target_fps = target_fps;

        print_banner();
        std::cout << "[WAKE] Launching scenario: \033[1;32m" << meta.name << "\033[0m\n"
                  << "[WAKE] Max steps: " << config.max_steps << " | Export directory: " << config.export_dir << "\n\n";

        auto runner = std::make_shared<wake::SimulationRunner>(std::move(scenario), config);
        g_active_runner = runner;

        std::signal(SIGINT, signal_handler);
        std::signal(SIGTERM, signal_handler);

        runner->run();

        std::cout << "\n[WAKE] Simulation complete. Average Throughput: "
                  << std::fixed << std::setprecision(2) << runner->average_mlups() << " MLUPS.\n"
                  << "[WAKE] Output artifacts saved to: " << config.export_dir << "\n";
        return 0;
    }

    if (cmd == "mesh") {
        if (argc < 3) {
            std::cerr << "Error: 'mesh' requires an input .obj file path.\n";
            return 1;
        }

        std::string obj_file = argv[2];
        int grid_size = 64;
        std::string output_vti = "voxelized_mesh.vti";

        for (int i = 3; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--grid" && i + 1 < argc) {
                grid_size = std::stoi(argv[++i]);
            } else if (arg == "--out" && i + 1 < argc) {
                output_vti = argv[++i];
            }
        }

        print_banner();
        std::cout << "[WAKE MESH] Loading 3D OBJ: " << obj_file << "\n";

        wake::ObjLoader loader;
        if (!loader.load_file(obj_file)) {
            std::cerr << "Error: Failed to open or parse OBJ file: " << obj_file << "\n";
            return 1;
        }

        std::cout << "[WAKE MESH] Loaded " << loader.vertex_count() << " vertices and "
                  << loader.triangle_count() << " triangles.\n";
        std::cout << "[WAKE MESH] Voxelizing into " << grid_size << "^3 structured lattice grid...\n";

        wake::LatticeGrid3D grid(grid_size, grid_size, grid_size, 0.02f);
        grid.initialize();

        wake::Vec3 pad(2.0f, 2.0f, 2.0f);
        wake::Vec3 target_max(grid_size - 3, grid_size - 3, grid_size - 3);
        loader.fit_into_box(pad, target_max);
        loader.voxelize_into_grid(grid, true);

        wake::VtkExporter::export_vti_3d(output_vti, grid);
        std::cout << "[WAKE MESH] Exported voxelized geometry to: " << output_vti << "\n";
        return 0;
    }

    std::cerr << "Error: Unknown command '" << cmd << "'. Use '" << argv[0] << " help' for usage.\n";
    return 1;
}
