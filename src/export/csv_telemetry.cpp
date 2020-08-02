#include "wake/export/csv_telemetry.hpp"
#include <iomanip>
#include <numeric>
#include <cmath>

namespace wake {

void TelemetryExporter::record_aero(const TelemetrySnapshot& sample) {
    aero_samples_.push_back(sample);
}

void TelemetryExporter::record_performance(const PerformanceMetric& perf) {
    perf_samples_.push_back(perf);
}

void TelemetryExporter::clear() {
    aero_samples_.clear();
    perf_samples_.clear();
}

bool TelemetryExporter::write_aero_csv(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        return false;
    }

    out << "step,time,fx,fy,fz,cd,cl,cm\n";
    out << std::fixed << std::setprecision(6);

    for (const auto& s : aero_samples_) {
        out << s.step << ","
            << s.time << ","
            << s.forces.fx << ","
            << s.forces.fy << ","
            << s.forces.fz << ","
            << s.forces.cd << ","
            << s.forces.cl << ","
            << s.forces.cm << "\n";
    }

    return true;
}

bool TelemetryExporter::write_performance_csv(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        return false;
    }

    out << "step,step_time_ms,mlups,max_speed,total_mass,kinetic_energy\n";
    out << std::fixed << std::setprecision(6);

    for (const auto& p : perf_samples_) {
        out << p.step << ","
            << p.step_time_ms << ","
            << p.mlups << ","
            << p.max_speed << ","
            << p.total_mass << ","
            << p.kinetic_energy << "\n";
    }

    return true;
}

bool TelemetryExporter::write_summary_json(
    const std::string& filepath,
    const std::string& simulation_name,
    real_t strouhal_number,
    real_t shedding_frequency
) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        return false;
    }

    // Statistical aggregation of aerodynamics
    double mean_cd = 0.0;
    double max_cd = -1e9;
    double min_cd = 1e9;
    double rms_cl = 0.0;

    if (!aero_samples_.empty()) {
        double sum_cd = 0.0;
        double sum_cl_sq = 0.0;
        for (const auto& s : aero_samples_) {
            sum_cd += s.forces.cd;
            sum_cl_sq += s.forces.cl * s.forces.cl;
            if (s.forces.cd > max_cd) max_cd = s.forces.cd;
            if (s.forces.cd < min_cd) min_cd = s.forces.cd;
        }
        mean_cd = sum_cd / aero_samples_.size();
        rms_cl = std::sqrt(sum_cl_sq / aero_samples_.size());
    } else {
        max_cd = 0.0;
        min_cd = 0.0;
    }

    // Average MLUPS
    double avg_mlups = 0.0;
    if (!perf_samples_.empty()) {
        double sum_mlups = 0.0;
        for (const auto& p : perf_samples_) {
            sum_mlups += p.mlups;
        }
        avg_mlups = sum_mlups / perf_samples_.size();
    }

    out << "{\n";
    out << "  \"simulation\": \"" << simulation_name << "\",\n";
    out << "  \"total_samples\": " << aero_samples_.size() << ",\n";
    out << "  \"performance\": {\n";
    out << "    \"average_mlups\": " << avg_mlups << "\n";
    out << "  },\n";
    out << "  \"aerodynamics\": {\n";
    out << "    \"mean_drag_coefficient\": " << mean_cd << ",\n";
    out << "    \"max_drag_coefficient\": " << max_cd << ",\n";
    out << "    \"min_drag_coefficient\": " << min_cd << ",\n";
    out << "    \"rms_lift_coefficient\": " << rms_cl << ",\n";
    out << "    \"shedding_frequency\": " << shedding_frequency << ",\n";
    out << "    \"strouhal_number\": " << strouhal_number << "\n";
    out << "  }\n";
    out << "}\n";

    return true;
}

} // namespace wake
