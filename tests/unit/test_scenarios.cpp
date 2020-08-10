#include "wake/scenarios/scenario_registry.hpp"
#include "wake/scenarios/karman_vortex.hpp"
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

int main() {
    test_scenario_registry();
    test_karman_vortex_execution();
    std::cout << "All scenario tests passed successfully!" << std::endl;
    return 0;
}
