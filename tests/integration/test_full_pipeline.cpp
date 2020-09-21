#include "wake/scenarios/scenario_registry.hpp"
#include "wake/runner/simulation_runner.hpp"
#include <iostream>
#include <cassert>
#include <filesystem>
#include <fstream>

using namespace wake;

int main() {
    std::cout << "[INTEGRATION TEST] Running full CFD simulation and export pipeline...\n";

    const std::string test_dir = "./integration_test_artifacts";
    std::filesystem::remove_all(test_dir);

    // 1. Instantiate scenario via registry
    auto scenario = ScenarioRegistry::instance().create("karman");
    assert(scenario != nullptr);

    // 2. Configure SimulationRunner
    RunnerConfig config;
    config.max_steps = 60;
    config.export_interval = 20;
    config.render_interval = 1000;
    config.live_terminal_render = false;
    config.export_dir = test_dir;
    config.target_fps = 0.0;

    SimulationRunner runner(std::move(scenario), config);
    runner.run();

    assert(runner.current_step() == 60);

    // 3. Verify exported artifacts
    // Check VTK VTI files exist
    assert(std::filesystem::exists(test_dir + "/karman_000000.vti"));
    assert(std::filesystem::exists(test_dir + "/karman_000001.vti"));

    // Check Netpbm PPM files exist
    assert(std::filesystem::exists(test_dir + "/karman_vorticity_000000.ppm"));
    assert(std::filesystem::exists(test_dir + "/karman_vorticity_000001.ppm"));

    // Check CSV telemetry exists and has content
    std::string aero_csv_path = test_dir + "/karman_aerodynamics.csv";
    assert(std::filesystem::exists(aero_csv_path));
    std::ifstream aero_file(aero_csv_path);
    std::string line;
    std::getline(aero_file, line);
    assert(line == "step,time,fx,fy,fz,cd,cl,cm");
    aero_file.close();

    // Check JSON summary exists and has valid structure
    std::string summary_json_path = test_dir + "/karman_summary.json";
    assert(std::filesystem::exists(summary_json_path));
    std::ifstream json_file(summary_json_path);
    std::string json_content((std::istreambuf_iterator<char>(json_file)), std::istreambuf_iterator<char>());
    assert(json_content.find("\"simulation\"") != std::string::npos);
    assert(json_content.find("\"aerodynamics\"") != std::string::npos);
    json_file.close();

    // 4. Clean up test artifacts
    std::filesystem::remove_all(test_dir);

    std::cout << "[INTEGRATION TEST] Full CFD pipeline PASSED!\n";
    return 0;
}
