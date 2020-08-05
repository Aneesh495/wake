#include "wake/scenarios/scenario_registry.hpp"
#include <stdexcept>

namespace wake {

ScenarioRegistry& ScenarioRegistry::instance() {
    static ScenarioRegistry inst;
    return inst;
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
