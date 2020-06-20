#include "wake/forces/aerodynamic_telemetry.hpp"
#include <fstream>
#include <cmath>
#include <numeric>
#include <iomanip>

namespace wake {

AerodynamicTelemetry::AerodynamicTelemetry(real_t characteristic_area, real_t chord)
    : characteristic_area_(characteristic_area), chord_(chord) {}

void AerodynamicTelemetry::record(uint64_t step, real_t time, const AerodynamicForces& forces) {
    history_.push_back({step, time, forces});
}

real_t AerodynamicTelemetry::compute_mean_cd(std::size_t start_step) const {
    if (history_.empty() || start_step >= history_.size()) return 0.0f;
    real_t sum = 0.0f;
    for (std::size_t i = start_step; i < history_.size(); ++i) {
        sum += history_[i].forces.cd;
    }
    return sum / static_cast<real_t>(history_.size() - start_step);
}

real_t AerodynamicTelemetry::compute_mean_cl(std::size_t start_step) const {
    if (history_.empty() || start_step >= history_.size()) return 0.0f;
    real_t sum = 0.0f;
    for (std::size_t i = start_step; i < history_.size(); ++i) {
        sum += history_[i].forces.cl;
    }
    return sum / static_cast<real_t>(history_.size() - start_step);
}

real_t AerodynamicTelemetry::compute_rms_cl(std::size_t start_step) const {
    if (history_.empty() || start_step >= history_.size()) return 0.0f;
    real_t mean = compute_mean_cl(start_step);
    real_t sum_sq = 0.0f;
    for (std::size_t i = start_step; i < history_.size(); ++i) {
        real_t diff = history_[i].forces.cl - mean;
        sum_sq += diff * diff;
    }
    return std::sqrt(sum_sq / static_cast<real_t>(history_.size() - start_step));
}

void AerodynamicTelemetry::export_csv(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) return;

    out << "step,time,fx,fy,fz,cd,cl,cm\n";
    out << std::fixed << std::setprecision(6);
    for (const auto& s : history_) {
        out << s.step << "," << s.time << ","
            << s.forces.fx << "," << s.forces.fy << "," << s.forces.fz << ","
            << s.forces.cd << "," << s.forces.cl << "," << s.forces.cm << "\n";
    }
}

void AerodynamicTelemetry::clear() {
    history_.clear();
}

} // namespace wake
