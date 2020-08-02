#pragma once

#include "wake/core/types.hpp"
#include "wake/forces/aerodynamic_telemetry.hpp"
#include <string>
#include <vector>
#include <fstream>

namespace wake {

struct PerformanceMetric {
    uint64_t step{0};
    double step_time_ms{0.0};
    double mlups{0.0};
    real_t max_speed{0.0f};
    real_t total_mass{0.0f};
    real_t kinetic_energy{0.0f};
};

class TelemetryExporter {
public:
    TelemetryExporter() = default;

    // Record aerodynamic sample
    void record_aero(const TelemetrySnapshot& sample);

    // Record performance benchmark sample
    void record_performance(const PerformanceMetric& perf);

    // Write full aerodynamic telemetry to CSV
    bool write_aero_csv(const std::string& filepath) const;

    // Write performance metrics to CSV
    bool write_performance_csv(const std::string& filepath) const;

    // Write combined simulation report to JSON format
    bool write_summary_json(
        const std::string& filepath,
        const std::string& simulation_name,
        real_t strouhal_number = 0.0f,
        real_t shedding_frequency = 0.0f
    ) const;

    // Reset recorded histories
    void clear();

    const std::vector<TelemetrySnapshot>& aero_samples() const { return aero_samples_; }
    const std::vector<PerformanceMetric>& perf_samples() const { return perf_samples_; }

private:
    std::vector<TelemetrySnapshot> aero_samples_;
    std::vector<PerformanceMetric> perf_samples_;
};

} // namespace wake
