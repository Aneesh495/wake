#include "wake/scenarios/scenario_registry.hpp"
#include "wake/scenarios/karman_vortex.hpp"
#include "wake/scenarios/naca_airfoil_stall.hpp"
#include "wake/scenarios/lid_driven_cavity.hpp"
#include "wake/scenarios/sphere_3d_wake.hpp"
#include "wake/scenarios/thermal_convection.hpp"
#include <stdexcept>

namespace wake {

ScenarioRegistry& ScenarioRegistry::instance() {
    static ScenarioRegistry inst;
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        initialize_builtin_scenarios();
    }
    return inst;
}

void ScenarioRegistry::initialize_builtin_scenarios() {
    auto& reg = ScenarioRegistry::instance();
    if (!reg.has_scenario("karman")) {
        reg.register_scenario("karman", []() -> std::unique_ptr<Scenario> {
            return std::make_unique<KarmanVortexScenario>();
        });
    }
    if (!reg.has_scenario("airfoil")) {
        reg.register_scenario("airfoil", []() -> std::unique_ptr<Scenario> {
            return std::make_unique<NacaAirfoilScenario>();
        });
    }
    if (!reg.has_scenario("cavity")) {
        reg.register_scenario("cavity", []() -> std::unique_ptr<Scenario> {
            return std::make_unique<LidDrivenCavityScenario>();
        });
    }
    if (!reg.has_scenario("sphere")) {
        reg.register_scenario("sphere", []() -> std::unique_ptr<Scenario> {
            return std::make_unique<Sphere3DWakeScenario>();
        });
    }
    if (!reg.has_scenario("thermal")) {
        reg.register_scenario("thermal", []() -> std::unique_ptr<Scenario> {
            return std::make_unique<ThermalConvectionScenario>();
        });
    }
}

void ScenarioRegistry::register_scenario(const std::string& name, ScenarioCreator creator) {
    registry_[name] = std::move(creator);
}

std::unique_ptr<Scenario> ScenarioRegistry::create(const std::string& name) const {
    auto it = registry_.find(name);
    if (it == registry_.end()) {
        throw std::runtime_error("Scenario not found in registry: " + name);
    }
    return (it->second)();
}

std::vector<std::string> ScenarioRegistry::list_scenarios() const {
    std::vector<std::string> names;
    names.reserve(registry_.size());
    for (const auto& [name, _] : registry_) {
        names.push_back(name);
    }
    return names;
}

bool ScenarioRegistry::has_scenario(const std::string& name) const {
    return registry_.find(name) != registry_.end();
}

} // namespace wake
