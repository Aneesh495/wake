#include "wake/core/types.hpp"
#include "wake/core/memory.hpp"
#include "wake/core/constants.hpp"
#include "wake/core/lattice_models.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>

namespace ispc {
extern "C" {
    void d2q9_compute_macroscopic_moments(
        const float f[],
        float rho[],
        float ux[],
        float uy[],
        int count
    );

    void d2q9_collide_bgk(
        const float f_in[],
        float f_out[],
        const uint8_t cell_type[],
        const float rho[],
        const float ux[],
        const float uy[],
        float omega,
        int count
    );
}
}

namespace wake {

// Scalar C++ reference implementation for D2Q9 moments
void scalar_d2q9_moments(
    const float f[],
    const uint8_t cell_type[],
    float rho[],
    float ux[],
    float uy[],
    float speed[],
    int count
) {
    for (int i = 0; i < count; ++i) {
        if (cell_type[i] == static_cast<uint8_t>(CellType::Solid)) {
            rho[i] = RHO_0;
            ux[i] = 0.0f;
            uy[i] = 0.0f;
            speed[i] = 0.0f;
            continue;
        }

        float r = 0.0f;
        float mx = 0.0f;
        float my = 0.0f;

        for (int q = 0; q < D2Q9::Q; ++q) {
            float fq = f[q * count + i];
            r += fq;
            mx += fq * static_cast<float>(D2Q9::CX[q]);
            my += fq * static_cast<float>(D2Q9::CY[q]);
        }

        float inv_rho = 1.0f / r;
        float u = mx * inv_rho;
        float v = my * inv_rho;

        rho[i] = r;
        ux[i] = u;
        uy[i] = v;
        speed[i] = std::sqrt(u * u + v * v);
    }
}

// Scalar C++ reference implementation for D2Q9 BGK collision
void scalar_d2q9_collide_bgk(
    const float f_in[],
    float f_out[],
    const uint8_t cell_type[],
    const float rho[],
    const float ux[],
    const float uy[],
    float omega,
    int count
) {
    for (int i = 0; i < count; ++i) {
        if (cell_type[i] == static_cast<uint8_t>(CellType::Solid)) {
            for (int q = 0; q < D2Q9::Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        float r = rho[i];
        float u = ux[i];
        float v = uy[i];
        float u_sq = u * u + v * v;

        for (int q = 0; q < D2Q9::Q; ++q) {
            float ci_u = static_cast<float>(D2Q9::CX[q]) * u + static_cast<float>(D2Q9::CY[q]) * v;
            float feq = D2Q9::W[q] * r * (1.0f + 3.0f * ci_u + 4.5f * ci_u * ci_u - 1.5f * u_sq);
            float f_curr = f_in[q * count + i];
            f_out[q * count + i] = f_curr - omega * (f_curr - feq);
        }
    }
}

} // namespace wake

int main() {
    using namespace wake;

    std::cout << "\n================================================================================\n";
    std::cout << "  WAKE CFD: D2Q9 ISPC SPMD vs Scalar C++ Performance Benchmark\n";
    std::cout << "================================================================================\n";
    std::cout << "  Architecture: arm64 / NEON (ISPC neon-i32x4)\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << std::setw(12) << "Grid Size"
              << std::setw(14) << "Total Nodes"
              << std::setw(18) << "Scalar (MLUPS)"
              << std::setw(18) << "ISPC (MLUPS)"
              << std::setw(14) << "Speedup"
              << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";

    const std::vector<int> resolutions = {128, 256, 512, 1024};
    const int iterations = 100;
    const float omega = 1.0f / 0.53f;

    for (int res : resolutions) {
        int count = res * res;
        AlignedBuffer<float> f_in(count * D2Q9::Q);
        AlignedBuffer<float> f_out_scalar(count * D2Q9::Q);
        AlignedBuffer<float> f_out_ispc(count * D2Q9::Q);
        AlignedBuffer<uint8_t> cell_type(count);

        AlignedBuffer<float> rho(count);
        AlignedBuffer<float> ux(count);
        AlignedBuffer<float> uy(count);
        AlignedBuffer<float> speed(count);

        for (int i = 0; i < count; ++i) {
            cell_type[i] = (i % 25 == 0) ? static_cast<uint8_t>(CellType::Solid) : static_cast<uint8_t>(CellType::Fluid);
            for (int q = 0; q < D2Q9::Q; ++q) {
                f_in[q * count + i] = D2Q9::W[q] * (1.0f + 0.02f * std::sin(static_cast<float>(i)));
            }
        }

        // 1. Benchmark Scalar C++
        auto t0_scalar = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iterations; ++it) {
            scalar_d2q9_moments(f_in.data(), cell_type.data(), rho.data(), ux.data(), uy.data(), speed.data(), count);
            scalar_d2q9_collide_bgk(f_in.data(), f_out_scalar.data(), cell_type.data(), rho.data(), ux.data(), uy.data(), omega, count);
        }
        auto t1_scalar = std::chrono::high_resolution_clock::now();
        double sec_scalar = std::chrono::duration<double>(t1_scalar - t0_scalar).count();
        double mlups_scalar = (static_cast<double>(count) * iterations / 1e6) / sec_scalar;

        // 2. Benchmark ISPC SPMD
        auto t0_ispc = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iterations; ++it) {
            ispc::d2q9_compute_macroscopic_moments(f_in.data(), rho.data(), ux.data(), uy.data(), count);
            ispc::d2q9_collide_bgk(f_in.data(), f_out_ispc.data(), cell_type.data(), rho.data(), ux.data(), uy.data(), omega, count);
        }
        auto t1_ispc = std::chrono::high_resolution_clock::now();
        double sec_ispc = std::chrono::duration<double>(t1_ispc - t0_ispc).count();
        double mlups_ispc = (static_cast<double>(count) * iterations / 1e6) / sec_ispc;

        double speedup = mlups_ispc / mlups_scalar;

        std::ostringstream res_str;
        res_str << res << "x" << res;

        std::cout << std::setw(12) << res_str.str()
                  << std::setw(14) << count
                  << std::setw(18) << std::fixed << std::setprecision(1) << mlups_scalar
                  << std::setw(18) << std::fixed << std::setprecision(1) << mlups_ispc
                  << std::setw(13) << std::fixed << std::setprecision(2) << speedup << "x"
                  << "\n";
    }

    std::cout << "================================================================================\n\n";
    return 0;
}
