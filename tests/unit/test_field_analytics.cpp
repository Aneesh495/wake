#include "wake/analytics/field_analytics.hpp"
#include <iostream>
#include <cmath>

using namespace wake;

int main() {
    std::cout << "[TEST] Running Flow Field Analytics (vorticity and Q-criterion) unit test...\n";

    // 1. Test 2D vorticity on rigid body rotation: u = (-omega * y, omega * x)
    // Analytical curl = 2 * omega
    const int NX = 16, NY = 16;
    const real_t omega_rot = 0.05f;
    const real_t expected_vorticity = 2.0f * omega_rot;

    LatticeGrid2D grid2d(NX, NY);
    grid2d.initialize_custom([&](int x, int y, real_t& rho, real_t& ux, real_t& uy) {
        rho = 1.0f;
        real_t rx = static_cast<real_t>(x) - 7.5f;
        real_t ry = static_cast<real_t>(y) - 7.5f;
        ux = -omega_rot * ry;
        uy = omega_rot * rx;
    });

    AlignedBuffer<real_t> vort_2d;
    FieldAnalytics::compute_vorticity(grid2d, vort_2d);

    // Check interior nodes (avoid periodic wrapping boundary discontinuity)
    for (int y = 2; y < NY - 2; ++y) {
        for (int x = 2; x < NX - 2; ++x) {
            real_t v = vort_2d[y * NX + x];
            if (std::abs(v - expected_vorticity) > 1e-5f) {
                std::cerr << "2D vorticity mismatch at (" << x << "," << y << "): expected "
                          << expected_vorticity << ", got " << v << "\n";
                return 1;
            }
        }
    }

    // 2. Test 3D Q-criterion on pure shear flow: ux = gamma * y, uy = 0, uz = 0
    // In pure shear, ||Omega||^2 = ||S||^2 -> Q = 0 exactly
    const int L = 16;
    LatticeGrid3D grid3d(L, L, L);
    const real_t shear_rate = 0.04f;

    grid3d.initialize_custom([&](int, int y, int, real_t& rho, real_t& ux, real_t& uy, real_t& uz) {
        rho = 1.0f;
        ux = shear_rate * static_cast<real_t>(y);
        uy = 0.0f;
        uz = 0.0f;
    });

    AlignedBuffer<real_t> wx, wy, wz, q_crit;
    FieldAnalytics::compute_3d_features(grid3d, wx, wy, wz, q_crit);

    for (int z = 2; z < L - 2; ++z) {
        for (int y = 2; y < L - 2; ++y) {
            for (int x = 2; x < L - 2; ++x) {
                int idx = (z * L + y) * L + x;
                // In pure simple shear, Q must be zero to machine precision
                if (std::abs(q_crit[idx]) > 1e-5f) {
                    std::cerr << "Pure shear Q-criterion non-zero: " << q_crit[idx] << "\n";
                    return 1;
                }
                // Vorticity wz should be -shear_rate
                if (std::abs(wz[idx] - (-shear_rate)) > 1e-5f) {
                    std::cerr << "Pure shear vorticity wz mismatch: expected "
                              << -shear_rate << ", got " << wz[idx] << "\n";
                    return 1;
                }
            }
        }
    }

    std::cout << "[PASS] Flow field analytics vorticity and Q-criterion verified against analytical tensor invariants.\n";
    return 0;
}
