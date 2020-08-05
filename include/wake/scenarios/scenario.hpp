#pragma once

#include "wake/core/types.hpp"
#include <string>
#include <memory>
#include <ostream>

namespace wake {

enum class SimulationDimension {
    Dim2D,
    Dim3D
};

struct ScenarioMetadata {
    std::string name;
    std::string description;
    SimulationDimension dimension{SimulationDimension::Dim2D};
    real_t reynolds_number{100.0f};
    real_t mach_number{0.05f};
    uint64_t default_steps{10000};
    uint64_t export_interval{100};
};

class Scenario {
public:
    virtual ~Scenario() = default;

    virtual const ScenarioMetadata& metadata() const noexcept = 0;

    // Simulation lifecycle
    virtual void setup() = 0;
    virtual void step() = 0;

    // State inspection
    virtual uint64_t current_step() const noexcept = 0;
    virtual real_t current_time() const noexcept = 0;
    virtual real_t current_mlups() const noexcept = 0;

    // Visualization and Exporting
    virtual void render_terminal(std::ostream& out) = 0;
    virtual void export_frame(const std::string& output_dir, int frame_idx) = 0;
    virtual void export_telemetry(const std::string& output_dir) = 0;

    virtual bool is_finished() const noexcept = 0;
};

} // namespace wake
