#include "wake/grid/lattice_grid_2d.hpp"
#include <iostream>
#include <cmath>
#include <vector>

namespace ispc {
extern "C" {
    void d2q9_collide_bgk(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float tau, float gx, float gy, int count);
    void d2q9_collide_trt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float tau_s, float tau_a, float gx, float gy, int count);
    void d2q9_collide_mrt(
        const float f_in[], float f_out[], const uint8_t cell_type[],
        float rho[], float ux[], float uy[], float tau_shear, float tau_bulk, float gx, float gy, int count);
}
}

int main() {
    std::cout << "[TEST] Running D2Q9 collision conservation (BGK, TRT, MRT) unit tests...\n";

    const int count = 32;
    wake::AlignedBuffer<float> f_in(count * wake::D2Q9::Q);
    wake::AlignedBuffer<float> f_out(count * wake::D2Q9::Q);
    wake::AlignedBuffer<uint8_t> cell_type(count);
    wake::AlignedBuffer<float> rho(count);
    wake::AlignedBuffer<float> ux(count);
    wake::AlignedBuffer<float> uy(count);

    // Node 0..30 fluid, node 31 solid
    for (int i = 0; i < count - 1; ++i) {
        cell_type[i] = static_cast<uint8_t>(wake::CellType::Fluid);
    }
    cell_type[count - 1] = static_cast<uint8_t>(wake::CellType::Solid);

    // Initial non-equilibrium state with perturbation
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        for (int i = 0; i < count; ++i) {
            float base = wake::D2Q9::W[q] * 1.0f;
            float pert = 0.01f * std::sin(static_cast<float>(q + i));
            f_in[q * count + i] = base + pert;
        }
    }

    // Compute pre-collision density and momentum for node 0
    float rho_pre = 0.0f;
    float jx_pre = 0.0f;
    float jy_pre = 0.0f;
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        float fval = f_in[q * count + 0];
        rho_pre += fval;
        jx_pre += fval * static_cast<float>(wake::D2Q9::CX[q]);
        jy_pre += fval * static_cast<float>(wake::D2Q9::CY[q]);
    }

    // 1. Test BGK collision conservation
    ispc::d2q9_collide_bgk(
        f_in.data(), f_out.data(), cell_type.data(),
        rho.data(), ux.data(), uy.data(), 0.7f, 0.0f, 0.0f, count
    );

    float rho_post_bgk = 0.0f;
    float jx_post_bgk = 0.0f;
    float jy_post_bgk = 0.0f;
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        float fval = f_out[q * count + 0];
        rho_post_bgk += fval;
        jx_post_bgk += fval * static_cast<float>(wake::D2Q9::CX[q]);
        jy_post_bgk += fval * static_cast<float>(wake::D2Q9::CY[q]);
    }

    if (std::abs(rho_post_bgk - rho_pre) > 1e-5f ||
        std::abs(jx_post_bgk - jx_pre) > 1e-5f ||
        std::abs(jy_post_bgk - jy_pre) > 1e-5f) {
        std::cerr << "BGK conservation failed: delta_rho=" << std::abs(rho_post_bgk - rho_pre)
                  << ", delta_jx=" << std::abs(jx_post_bgk - jx_pre) << "\n";
        return 1;
    }

    // 2. Test TRT collision conservation
    ispc::d2q9_collide_trt(
        f_in.data(), f_out.data(), cell_type.data(),
        rho.data(), ux.data(), uy.data(), 0.7f, 1.2f, 0.0f, 0.0f, count
    );

    float rho_post_trt = 0.0f;
    float jx_post_trt = 0.0f;
    float jy_post_trt = 0.0f;
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        float fval = f_out[q * count + 0];
        rho_post_trt += fval;
        jx_post_trt += fval * static_cast<float>(wake::D2Q9::CX[q]);
        jy_post_trt += fval * static_cast<float>(wake::D2Q9::CY[q]);
    }

    if (std::abs(rho_post_trt - rho_pre) > 1e-5f ||
        std::abs(jx_post_trt - jx_pre) > 1e-5f ||
        std::abs(jy_post_trt - jy_pre) > 1e-5f) {
        std::cerr << "TRT conservation failed: delta_rho=" << std::abs(rho_post_trt - rho_pre) << "\n";
        return 1;
    }

    // 3. Test MRT collision conservation
    ispc::d2q9_collide_mrt(
        f_in.data(), f_out.data(), cell_type.data(),
        rho.data(), ux.data(), uy.data(), 0.7f, 1.6f, 0.0f, 0.0f, count
    );

    float rho_post_mrt = 0.0f;
    float jx_post_mrt = 0.0f;
    float jy_post_mrt = 0.0f;
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        float fval = f_out[q * count + 0];
        rho_post_mrt += fval;
        jx_post_mrt += fval * static_cast<float>(wake::D2Q9::CX[q]);
        jy_post_mrt += fval * static_cast<float>(wake::D2Q9::CY[q]);
    }

    if (std::abs(rho_post_mrt - rho_pre) > 1e-5f ||
        std::abs(jx_post_mrt - jx_pre) > 1e-5f ||
        std::abs(jy_post_mrt - jy_pre) > 1e-5f) {
        std::cerr << "MRT conservation failed: delta_rho=" << std::abs(rho_post_mrt - rho_pre) << "\n";
        return 1;
    }

    // 4. Verify solid node invariance
    for (int q = 0; q < wake::D2Q9::Q; ++q) {
        float diff = std::abs(f_out[q * count + (count - 1)] - f_in[q * count + (count - 1)]);
        if (diff > 1e-6f) {
            std::cerr << "Solid node altered by collision operator! diff=" << diff << "\n";
            return 1;
        }
    }

    std::cout << "[PASS] D2Q9 collision conservation tests passed for BGK, TRT, and MRT.\n";
    return 0;
}
