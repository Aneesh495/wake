#pragma once

#include "wake/scenarios/scenario.hpp"
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_map>

namespace wake {

using ScenarioCreator = std::function<std::unique_ptr<Scenario>()>;

class ScenarioRegistry {
public:
    static ScenarioRegistry& instance();

    // Register a new simulation scenario factory
    void register_scenario(const std::string& name, ScenarioCreator creator);

    // Create scenario by name
    std::unique_ptr<Scenario> create(const std::string& name) const;

    // List all registered scenario names
    std::vector<std::string> list_scenarios() const;

    // Check if scenario is registered
    bool has_scenario(const std::string& name) const;

private:
    ScenarioRegistry() = default;
    std::unordered_map<std::string, ScenarioCreator> registry_;
};

// Helper macro for auto-registering scenario types
#define REGISTER_SCENARIO(NameStr, ScenarioClass) \
    namespace { \
        struct ScenarioClass##_Registrar { \
            ScenarioClass##_Registrar() { \
                ::wake::ScenarioRegistry::instance().register_scenario( \
                    NameStr, \
                    []() -> std::unique_ptr<::wake::Scenario> { \
                        return std::make_unique<ScenarioClass>(); \
                    } \
                ); \
            } \
        }; \
        static ScenarioClass##_Registrar s_##ScenarioClass##_registrar; \
    }

} // namespace wake
