#pragma once

#include "wake/forces/momentum_exchange.hpp"
#include <vector>
#include <string>

namespace wake {

// Telemetry record for a single simulation timestep
struct TelemetrySnapshot {
    uint64_t step;
    real_t time;
    AerodynamicForces forces;
};

// Records and aggregates time-series aerodynamic telemetry
class AerodynamicTelemetry {
public:
    explicit AerodynamicTelemetry(real_t characteristic_area = 1.0f, real_t chord = 1.0f);

    void record(uint64_t step, real_t time, const AerodynamicForces& forces);

    // Statistical queries across a specified step window
    real_t compute_mean_cd(std::size_t start_step = 0) const;
    real_t compute_mean_cl(std::size_t start_step = 0) const;
    real_t compute_rms_cl(std::size_t start_step = 0) const;

    real_t characteristic_area() const noexcept { return characteristic_area_; }
    real_t chord() const noexcept { return chord_; }

    const std::vector<TelemetrySnapshot>& history() const noexcept { return history_; }
    std::size_t size() const noexcept { return history_.size(); }

    void export_csv(const std::string& filepath) const;
    void clear();

private:
    real_t characteristic_area_{1.0f};
    real_t chord_{1.0f};
    std::vector<TelemetrySnapshot> history_;
};

} // namespace wake
