#include "wake/forces/momentum_exchange.hpp"
#include "wake/forces/strouhal_analyzer.hpp"
#include "wake/forces/aerodynamic_telemetry.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include <iostream>
#include <cmath>
#include <vector>

using namespace wake;

namespace ispc {
extern "C" {
    void d2q9_compute_momentum_exchange_forces(
        const float f_post_collision[],
        const uint8_t cell_type[],
        int nx, int ny,
        float center_x, float center_y,
        float &fx_out, float &fy_out, float &torque_z_out
    );
}
}

int main() {
    std::cout << "[TEST] Running Momentum Exchange Method (MEM) and FFT Strouhal analyzer unit test...\n";

    // 1. Test 2D Momentum Exchange Method kernel
    const int NX = 16, NY = 16;
    const int count = NX * NY;
    AlignedBuffer<float> f(count * D2Q9::Q);
    AlignedBuffer<uint8_t> cell_type(count);

    for (int i = 0; i < count; ++i) {
        cell_type[i] = static_cast<uint8_t>(CellType::Fluid);
        for (int q = 0; q < D2Q9::Q; ++q) {
            f[q * count + i] = D2Q9::W[q] * 1.0f;
        }
    }

    // Set a solid 2x2 square in the center: (x, y) in [7, 8] x [7, 8]
    for (int y = 7; y <= 8; ++y) {
        for (int x = 7; x <= 8; ++x) {
            cell_type[y * NX + x] = static_cast<uint8_t>(CellType::Solid);
        }
    }

    float fx = 0.0f, fy = 0.0f, tz = 0.0f;
    ispc::d2q9_compute_momentum_exchange_forces(
        f.data(), cell_type.data(), NX, NY, 7.5f, 7.5f, fx, fy, tz
    );

    // In static equilibrium with uniform distributions, net symmetric forces should be zero
    if (std::abs(fx) > 1e-5f || std::abs(fy) > 1e-5f) {
        std::cerr << "Equilibrium force balance failed: fx=" << fx << ", fy=" << fy << "\n";
        return 1;
    }

    // Introduce asymmetric momentum on west face of obstacle: node (6, 7) pointing east (q=1)
    f[1 * count + (7 * NX + 6)] += 0.05f;
    ispc::d2q9_compute_momentum_exchange_forces(
        f.data(), cell_type.data(), NX, NY, 7.5f, 7.5f, fx, fy, tz
    );

    // d_fx = 2.0 * 0.05 * 1.0 = 0.1
    if (std::abs(fx - 0.1f) > 1e-5f) {
        std::cerr << "Asymmetric force impulse failed: expected fx=0.1, got " << fx << "\n";
        return 1;
    }

    // 2. Test Strouhal FFT Frequency Analyzer on synthetic sinusoidal vortex signal
    StrouhalAnalyzer analyzer;
    const real_t true_frequency = 0.025f; // Hz (period T = 40 steps)
    const real_t dt = 1.0f;
    const std::size_t num_samples = 512;

    std::vector<real_t> signal(num_samples);
    for (std::size_t i = 0; i < num_samples; ++i) {
        real_t t = static_cast<real_t>(i) * dt;
        signal[i] = 0.35f * std::sin(TWO_PI * true_frequency * t) + 0.02f * std::cos(TWO_PI * 0.1f * t);
    }

    real_t detected_freq = analyzer.analyze_frequency(signal, dt);
    real_t freq_err = std::abs(detected_freq - true_frequency);

    std::cout << "Target Frequency: " << true_frequency << " Hz, Detected Frequency: " << detected_freq
              << " Hz (Error: " << freq_err << " Hz)\n";

    if (freq_err > 0.003f) {
        std::cerr << "FFT frequency extraction failed! Detected: " << detected_freq << "\n";
        return 1;
    }

    real_t cylinder_diameter = 20.0f;
    real_t u_inf = 0.1f;
    real_t strouhal = analyzer.compute_strouhal_number(detected_freq, cylinder_diameter, u_inf);
    real_t expected_st = (true_frequency * cylinder_diameter) / u_inf; // 0.025 * 20 / 0.1 = 5.0
    if (std::abs(strouhal - expected_st) > 0.6f) {
        std::cerr << "Strouhal calculation failed: " << strouhal << "\n";
        return 1;
    }

    std::cout << "[PASS] Momentum Exchange Method and Strouhal FFT analyzer passed successfully.\n";
    return 0;
}
