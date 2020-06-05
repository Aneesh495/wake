#include "wake/grid/lattice_grid_3d.hpp"
#include <iostream>
#include <cmath>

namespace ispc {
extern "C" {
    void d3q19_collide_bgk(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float uz[],
        float tau, float gx, float gy, float gz, int count);
    void d3q19_collide_trt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float uz[],
        float tau_s, float tau_a, float gx, float gy, float gz, int count);
    void d3q19_collide_mrt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float uz[],
        float tau_shear, float tau_bulk, float gx, float gy, float gz, int count);
}
}

int main() {
    std::cout << "[TEST] Running 3D D3Q19 collision conservation (BGK, TRT, MRT) unit tests...\n";

    const int count = 16;
    wake::AlignedBuffer<float> f_in(count * wake::D3Q19::Q);
    wake::AlignedBuffer<float> f_out(count * wake::D3Q19::Q);
    wake::AlignedBuffer<uint8_t> cell_type(count);
    wake::AlignedBuffer<float> rho(count);
    wake::AlignedBuffer<float> ux(count);
    wake::AlignedBuffer<float> uy(count);
    wake::AlignedBuffer<float> uz(count);

    for (int i = 0; i < count - 1; ++i) {
        cell_type[i] = static_cast<uint8_t>(wake::CellType::Fluid);
    }
    cell_type[count - 1] = static_cast<uint8_t>(wake::CellType::Solid);

    for (int q = 0; q < wake::D3Q19::Q; ++q) {
        for (int i = 0; i < count; ++i) {
            float base = wake::D3Q19::W[q] * 1.0f;
            float pert = 0.005f * std::sin(static_cast<float>(q * 3 + i));
            f_in[q * count + i] = base + pert;
        }
    }

    float rho_pre = 0.0f, jx_pre = 0.0f, jy_pre = 0.0f, jz_pre = 0.0f;
    for (int q = 0; q < wake::D3Q19::Q; ++q) {
        float fval = f_in[q * count + 0];
        rho_pre += fval;
        jx_pre += fval * static_cast<float>(wake::D3Q19::CX[q]);
        jy_pre += fval * static_cast<float>(wake::D3Q19::CY[q]);
        jz_pre += fval * static_cast<float>(wake::D3Q19::CZ[q]);
    }

    // 1. BGK
    ispc::d3q19_collide_bgk(
        f_in.data(), f_out.data(), cell_type.data(),
        rho.data(), ux.data(), uy.data(), uz.data(),
        0.75f, 0.0f, 0.0f, 0.0f, count
    );

    float rho_bgk = 0.0f, jx_bgk = 0.0f, jy_bgk = 0.0f, jz_bgk = 0.0f;
    for (int q = 0; q < wake::D3Q19::Q; ++q) {
        float fval = f_out[q * count + 0];
        rho_bgk += fval;
        jx_bgk += fval * static_cast<float>(wake::D3Q19::CX[q]);
        jy_bgk += fval * static_cast<float>(wake::D3Q19::CY[q]);
        jz_bgk += fval * static_cast<float>(wake::D3Q19::CZ[q]);
    }

    if (std::abs(rho_bgk - rho_pre) > 1e-5f ||
        std::abs(jx_bgk - jx_pre) > 1e-5f ||
        std::abs(jy_bgk - jy_pre) > 1e-5f ||
        std::abs(jz_bgk - jz_pre) > 1e-5f) {
        std::cerr << "D3Q19 BGK conservation failed!\n";
        return 1;
    }

    // 2. TRT
    ispc::d3q19_collide_trt(
        f_in.data(), f_out.data(), cell_type.data(),
        rho.data(), ux.data(), uy.data(), uz.data(),
        0.75f, 1.3f, 0.0f, 0.0f, 0.0f, count
    );

    float rho_trt = 0.0f, jx_trt = 0.0f;
    for (int q = 0; q < wake::D3Q19::Q; ++q) {
        float fval = f_out[q * count + 0];
        rho_trt += fval;
        jx_trt += fval * static_cast<float>(wake::D3Q19::CX[q]);
    }

    if (std::abs(rho_trt - rho_pre) > 1e-5f || std::abs(jx_trt - jx_pre) > 1e-5f) {
        std::cerr << "D3Q19 TRT conservation failed!\n";
        return 1;
    }

    // 3. MRT
    ispc::d3q19_collide_mrt(
        f_in.data(), f_out.data(), cell_type.data(),
        rho.data(), ux.data(), uy.data(), uz.data(),
        0.75f, 1.2f, 0.0f, 0.0f, 0.0f, count
    );

    float rho_mrt = 0.0f, jx_mrt = 0.0f;
    for (int q = 0; q < wake::D3Q19::Q; ++q) {
        float fval = f_out[q * count + 0];
        rho_mrt += fval;
        jx_mrt += fval * static_cast<float>(wake::D3Q19::CX[q]);
    }

    if (std::abs(rho_mrt - rho_pre) > 1e-5f || std::abs(jx_mrt - jx_pre) > 1e-5f) {
        std::cerr << "D3Q19 MRT conservation failed!\n";
        return 1;
    }

    std::cout << "[PASS] 3D D3Q19 collision conservation verified across BGK, TRT, and MRT.\n";
    return 0;
}
