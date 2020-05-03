#include "wake/grid/lattice_grid_2d.hpp"
#include <iostream>
#include <cmath>
#include <cassert>

namespace ispc {
extern "C" {
    void d2q9_compute_macroscopic_moments(
        const float f[], float rho[], float ux[], float uy[], int count);
    void d2q9_compute_velocity_magnitude(
        const float ux[], const float uy[], float speed[], int count);
}
}

int main() {
    std::cout << "[TEST] Running D2Q9 moments and equilibrium unit test...\n";

    const int NX = 16;
    const int NY = 16;
    const int count = NX * NY;

    wake::AlignedBuffer<float> f(count * wake::D2Q9::Q);
    wake::AlignedBuffer<float> rho(count);
    wake::AlignedBuffer<float> ux(count);
    wake::AlignedBuffer<float> uy(count);
    wake::AlignedBuffer<float> speed(count);

    const float target_rho = 1.05f;
    const float target_ux = 0.08f;
    const float target_uy = -0.04f;
    const float u_sq = target_ux * target_ux + target_uy * target_uy;

    // Populate f with exact analytical Maxwell-Boltzmann equilibrium
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        float ci_u = static_cast<float>(wake::D2Q9::CX[q]) * target_ux +
                     static_cast<float>(wake::D2Q9::CY[q]) * target_uy;
        float feq = wake::D2Q9::W[q] * target_rho * (
            1.0f + wake::INV_CS2 * ci_u + wake::INV_2_CS4 * ci_u * ci_u - wake::INV_2_CS2 * u_sq
        );
        for (int i = 0; i < count; ++i) {
            f[q * count + i] = feq;
        }
    }

    // Call ISPC kernel
    ispc::d2q9_compute_macroscopic_moments(
        f.data(), rho.data(), ux.data(), uy.data(), count
    );
    ispc::d2q9_compute_velocity_magnitude(
        ux.data(), uy.data(), speed.data(), count
    );

    // Verify precision across all nodes
    for (int i = 0; i < count; ++i) {
        float err_rho = std::abs(rho[i] - target_rho);
        float err_ux  = std::abs(ux[i] - target_ux);
        float err_uy  = std::abs(uy[i] - target_uy);
        float expected_speed = std::sqrt(u_sq);
        float err_speed = std::abs(speed[i] - expected_speed);

        if (err_rho > 1e-5f || err_ux > 1e-5f || err_uy > 1e-5f || err_speed > 1e-5f) {
            std::cerr << "FAILED at node " << i
                      << ": rho=" << rho[i] << " (err=" << err_rho << ")"
                      << ", ux=" << ux[i] << " (err=" << err_ux << ")"
                      << ", uy=" << uy[i] << " (err=" << err_uy << ")\n";
            return 1;
        }
    }

    std::cout << "[PASS] D2Q9 moments and equilibrium test passed successfully.\n";
    return 0;
}
