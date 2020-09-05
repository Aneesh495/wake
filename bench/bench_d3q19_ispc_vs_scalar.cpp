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
    void d3q19_compute_macroscopic_moments(
        const float f[],
        float rho[],
        float ux[],
        float uy[],
        float uz[],
        int count
    );

    void d3q19_collide_bgk(
        const float f_in[],
        float f_out[],
        const uint8_t cell_type[],
        const float rho[],
        const float ux[],
        const float uy[],
        const float uz[],
        float omega,
        float gx,
        float gy,
        float gz,
        int count
    );
}
}

namespace wake {

// Scalar C++ reference implementation for D3Q19 moments
void scalar_d3q19_moments(
    const float f[],
    const uint8_t cell_type[],
    float rho[],
    float ux[],
    float uy[],
    float uz[],
    float speed[],
    int count
) {
    for (int i = 0; i < count; ++i) {
        if (cell_type[i] == static_cast<uint8_t>(CellType::Solid)) {
            rho[i] = RHO_0;
            ux[i] = 0.0f;
            uy[i] = 0.0f;
            uz[i] = 0.0f;
            speed[i] = 0.0f;
            continue;
        }

        float r = 0.0f;
        float mx = 0.0f;
        float my = 0.0f;
        float mz = 0.0f;

        for (int q = 0; q < D3Q19::Q; ++q) {
            float fq = f[q * count + i];
            r += fq;
            mx += fq * static_cast<float>(D3Q19::CX[q]);
            my += fq * static_cast<float>(D3Q19::CY[q]);
            mz += fq * static_cast<float>(D3Q19::CZ[q]);
        }

        float inv_rho = 1.0f / r;
        float u = mx * inv_rho;
        float v = my * inv_rho;
        float w = mz * inv_rho;

        rho[i] = r;
        ux[i] = u;
        uy[i] = v;
        uz[i] = w;
        speed[i] = std::sqrt(u * u + v * v + w * w);
    }
}

// Scalar C++ reference implementation for D3Q19 BGK collision
void scalar_d3q19_collide_bgk(
    const float f_in[],
    float f_out[],
    const uint8_t cell_type[],
    const float rho[],
    const float ux[],
    const float uy[],
    const float uz[],
    float omega,
    float gx,
    float gy,
    float gz,
    int count
) {
    for (int i = 0; i < count; ++i) {
        if (cell_type[i] == static_cast<uint8_t>(CellType::Solid)) {
            for (int q = 0; q < D3Q19::Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        float r = rho[i];
        float u = ux[i];
        float v = uy[i];
        float w = uz[i];
        float u_sq = u * u + v * v + w * w;

        for (int q = 0; q < D3Q19::Q; ++q) {
            float cx = static_cast<float>(D3Q19::CX[q]);
            float cy = static_cast<float>(D3Q19::CY[q]);
            float cz = static_cast<float>(D3Q19::CZ[q]);
            float ci_u = cx * u + cy * v + cz * w;

            float feq = D3Q19::W[q] * r * (1.0f + 3.0f * ci_u + 4.5f * ci_u * ci_u - 1.5f * u_sq);
            float ci_g = cx * gx + cy * gy + cz * gz;
            float force_term = (1.0f - 0.5f * omega) * D3Q19::W[q] * (3.0f * ci_g + 9.0f * ci_u * ci_g - 3.0f * (u * gx + v * gy + w * gz));

            float f_curr = f_in[q * count + i];
            f_out[q * count + i] = f_curr - omega * (f_curr - feq) + force_term;
        }
    }
}

} // namespace wake

int main() {
    using namespace wake;

    std::cout << "\n================================================================================\n";
    std::cout << "  WAKE CFD: 3D D3Q19 ISPC SPMD vs Scalar C++ Performance Benchmark\n";
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

    const std::vector<int> resolutions = {32, 48, 64, 80};
    const int iterations = 50;
    const float omega = 1.0f / 0.53f;

    for (int res : resolutions) {
        int count = res * res * res;
        AlignedBuffer<float> f_in(count * D3Q19::Q);
        AlignedBuffer<float> f_out_scalar(count * D3Q19::Q);
        AlignedBuffer<float> f_out_ispc(count * D3Q19::Q);
        AlignedBuffer<uint8_t> cell_type(count);

        AlignedBuffer<float> rho(count);
        AlignedBuffer<float> ux(count);
        AlignedBuffer<float> uy(count);
        AlignedBuffer<float> uz(count);
        AlignedBuffer<float> speed(count);

        for (int i = 0; i < count; ++i) {
            cell_type[i] = (i % 30 == 0) ? static_cast<uint8_t>(CellType::Solid) : static_cast<uint8_t>(CellType::Fluid);
            for (int q = 0; q < D3Q19::Q; ++q) {
                f_in[q * count + i] = D3Q19::W[q] * (1.0f + 0.01f * std::sin(static_cast<float>(i)));
            }
        }

        // 1. Benchmark Scalar C++
        auto t0_scalar = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iterations; ++it) {
            scalar_d3q19_moments(f_in.data(), cell_type.data(), rho.data(), ux.data(), uy.data(), uz.data(), speed.data(), count);
            scalar_d3q19_collide_bgk(f_in.data(), f_out_scalar.data(), cell_type.data(), rho.data(), ux.data(), uy.data(), uz.data(), omega, 0.0f, 0.0f, 0.0f, count);
        }
        auto t1_scalar = std::chrono::high_resolution_clock::now();
        double sec_scalar = std::chrono::duration<double>(t1_scalar - t0_scalar).count();
        double mlups_scalar = (static_cast<double>(count) * iterations / 1e6) / sec_scalar;

        // 2. Benchmark ISPC SPMD
        auto t0_ispc = std::chrono::high_resolution_clock::now();
        for (int it = 0; it < iterations; ++it) {
            ispc::d3q19_compute_macroscopic_moments(f_in.data(), rho.data(), ux.data(), uy.data(), uz.data(), count);
            ispc::d3q19_collide_bgk(f_in.data(), f_out_ispc.data(), cell_type.data(), rho.data(), ux.data(), uy.data(), uz.data(), omega, 0.0f, 0.0f, 0.0f, count);
        }
        auto t1_ispc = std::chrono::high_resolution_clock::now();
        double sec_ispc = std::chrono::duration<double>(t1_ispc - t0_ispc).count();
        double mlups_ispc = (static_cast<double>(count) * iterations / 1e6) / sec_ispc;

        double speedup = mlups_ispc / mlups_scalar;

        std::ostringstream res_str;
        res_str << res << "^3";

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
