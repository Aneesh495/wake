#include "wake/scenarios/scenario_registry.hpp"
#include "wake/scenarios/karman_vortex.hpp"
#include "wake/scenarios/naca_airfoil_stall.hpp"
#include "wake/scenarios/lid_driven_cavity.hpp"
#include <iostream>
#include <cassert>

using namespace wake;

void test_scenario_registry() {
    std::cout << "Testing Scenario Registry..." << std::endl;

    auto list = ScenarioRegistry::instance().list_scenarios();
    assert(!list.empty());
    assert(ScenarioRegistry::instance().has_scenario("karman"));

    auto scenario = ScenarioRegistry::instance().create("karman");
    assert(scenario != nullptr);
    assert(scenario->metadata().name == "karman_vortex_street");

    std::cout << "Scenario Registry PASSED." << std::endl;
}

void test_karman_vortex_execution() {
    std::cout << "Testing Karman Vortex scenario execution (small grid)..." << std::endl;

    // Run small grid for 50 steps
    KarmanVortexScenario scenario(80, 40, 100.0f, 0.05f, 10.0f);
    scenario.setup();

    for (int step = 0; step < 50; ++step) {
        scenario.step();
    }

    assert(scenario.current_step() == 50);
    assert(scenario.current_mlups() > 0.0);
    // Drag coefficient should be positive and non-zero
    assert(scenario.latest_cd() > 0.0f);

    std::cout << "Karman Vortex execution PASSED (Cd = " << scenario.latest_cd()
              << ", MLUPS = " << scenario.current_mlups() << ")." << std::endl;
}

void test_airfoil_scenario_execution() {
    std::cout << "Testing NACA Airfoil scenario execution (small grid)..." << std::endl;

    NacaAirfoilScenario scenario("0012", 8.0f, 100, 50, 500.0f, 0.05f, 30.0f);
    scenario.setup();

    for (int step = 0; step < 50; ++step) {
        scenario.step();
    }

    assert(scenario.current_step() == 50);
    assert(scenario.current_mlups() > 0.0);
    assert(scenario.latest_cd() > 0.0f);

    std::cout << "NACA Airfoil execution PASSED (Cd = " << scenario.latest_cd()
              << ", Cl = " << scenario.latest_cl()
              << ", MLUPS = " << scenario.current_mlups() << ")." << std::endl;
}

void test_cavity_scenario_execution() {
    std::cout << "Testing 3D Lid-Driven Cavity scenario execution..." << std::endl;

    LidDrivenCavityScenario scenario(24, 100.0f, 0.08f);
    scenario.setup();

    for (int step = 0; step < 20; ++step) {
        scenario.step();
    }

    assert(scenario.current_step() == 20);
    assert(scenario.current_mlups() > 0.0);

    auto ux_prof = scenario.extract_vertical_centerline_ux();
    assert(ux_prof.size() == 24);

    std::cout << "3D Lid-Driven Cavity execution PASSED (MLUPS = " << scenario.current_mlups() << ")." << std::endl;
}

int main() {
    test_scenario_registry();
    test_karman_vortex_execution();
    test_airfoil_scenario_execution();
    test_cavity_scenario_execution();
    std::cout << "All scenario tests passed successfully!" << std::endl;
    return 0;
}
