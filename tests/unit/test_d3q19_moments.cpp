#include "wake/grid/lattice_grid_3d.hpp"
#include <iostream>
#include <cmath>

namespace ispc {
extern "C" {
    void d3q19_compute_macroscopic_moments(
        const float f[], float rho[], float ux[], float uy[], float uz[], int count);
    void d3q19_compute_velocity_magnitude(
        const float ux[], const float uy[], const float uz[], float speed[], int count);
}
}

int main() {
    std::cout << "[TEST] Running 3D D3Q19 moments and equilibrium unit test...\n";

    const int NX = 8, NY = 8, NZ = 8;
    const int count = NX * NY * NZ;

    wake::AlignedBuffer<float> f(count * wake::D3Q19::Q);
    wake::AlignedBuffer<float> rho(count);
    wake::AlignedBuffer<float> ux(count);
    wake::AlignedBuffer<float> uy(count);
    wake::AlignedBuffer<float> uz(count);
    wake::AlignedBuffer<float> speed(count);

    const float target_rho = 1.02f;
    const float target_ux = 0.05f;
    const float target_uy = -0.03f;
    const float target_uz = 0.04f;
    const float u_sq = target_ux * target_ux + target_uy * target_uy + target_uz * target_uz;

    for (int q = 0; q < wake::D3Q19::Q; ++q) {
        float ci_u = static_cast<float>(wake::D3Q19::CX[q]) * target_ux +
                     static_cast<float>(wake::D3Q19::CY[q]) * target_uy +
                     static_cast<float>(wake::D3Q19::CZ[q]) * target_uz;
        float feq = wake::D3Q19::W[q] * target_rho * (
            1.0f + wake::INV_CS2 * ci_u + wake::INV_2_CS4 * ci_u * ci_u - wake::INV_2_CS2 * u_sq
        );
        for (int i = 0; i < count; ++i) {
            f[q * count + i] = feq;
        }
    }

    ispc::d3q19_compute_macroscopic_moments(
        f.data(), rho.data(), ux.data(), uy.data(), uz.data(), count
    );
    ispc::d3q19_compute_velocity_magnitude(
        ux.data(), uy.data(), uz.data(), speed.data(), count
    );

    for (int i = 0; i < count; ++i) {
        float err_rho = std::abs(rho[i] - target_rho);
        float err_ux  = std::abs(ux[i] - target_ux);
        float err_uy  = std::abs(uy[i] - target_uy);
        float err_uz  = std::abs(uz[i] - target_uz);
        float exp_speed = std::sqrt(u_sq);
        float err_spd = std::abs(speed[i] - exp_speed);

        if (err_rho > 1e-5f || err_ux > 1e-5f || err_uy > 1e-5f || err_uz > 1e-5f || err_spd > 1e-5f) {
            std::cerr << "D3Q19 moment check failed at node " << i << "\n";
            return 1;
        }
    }

    std::cout << "[PASS] 3D D3Q19 moments and equilibrium verified successfully.\n";
    return 0;
}
