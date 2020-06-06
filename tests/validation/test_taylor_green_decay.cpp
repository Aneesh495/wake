#include "wake/grid/lattice_grid_3d.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace wake;

int main() {
    std::cout << "[VALIDATION] Running 3D Taylor-Green vortex kinetic energy decay validation test...\n";

    const int L = 32;
    const real_t tau = 0.8f;
    const real_t nu = CS2 * (tau - 0.5f); // 0.1
    const real_t u0 = 0.04f;
    const real_t k = TWO_PI / static_cast<real_t>(L);

    LatticeGrid3D grid(L, L, L, tau);
    grid.set_periodic(true);

    // Initialize with analytical Taylor-Green velocity and pressure fields
    grid.initialize_custom([&](int x, int y, int z, real_t& rho, real_t& ux, real_t& uy, real_t& uz) {
        real_t kx = k * static_cast<real_t>(x);
        real_t ky = k * static_cast<real_t>(y);
        real_t kz = k * static_cast<real_t>(z);

        ux = u0 * std::sin(kx) * std::cos(ky) * std::cos(kz);
        uy = -u0 * std::cos(kx) * std::sin(ky) * std::cos(kz);
        uz = 0.0f;

        real_t p_pert = (RHO_0 * u0 * u0 / 16.0f) * (std::cos(2.0f * kx) + std::cos(2.0f * ky)) * (std::cos(2.0f * kz) + 2.0f);
        rho = RHO_0 + p_pert / CS2;
    });

    real_t E0 = grid.compute_total_kinetic_energy();
    std::cout << "Initial Total Kinetic Energy E0: " << E0 << "\n";

    const int total_steps = 100;
    const int sample_interval = 20;

    std::cout << "  Step | Simulated E(t)/E0 | Analytical E(t)/E0 | Relative Diff\n";
    std::cout << "-------+-------------------+--------------------+---------------\n";

    for (int step = 1; step <= total_steps; ++step) {
        grid.step(CollisionType::MRT);

        if (step % sample_interval == 0) {
            real_t E_sim = grid.compute_total_kinetic_energy() / E0;
            // Analytical 3D laminar kinetic energy decay: exp(-2 * (kx^2 + ky^2 + kz^2) * nu * t) = exp(-6 * k^2 * nu * t)
            real_t t_sim = static_cast<real_t>(step);
            real_t E_ana = std::exp(-6.0f * k * k * nu * t_sim);

            real_t diff = std::abs(E_sim - E_ana);
            std::cout << std::setw(6) << step << " | "
                      << std::setw(17) << E_sim << " | "
                      << std::setw(18) << E_ana << " | "
                      << std::setw(13) << diff << "\n";

            if (diff > 0.025f) { // Within 2.5%
                std::cerr << "Taylor-Green decay test failed at step " << step << " with diff " << diff << "\n";
                return 1;
            }
        }
    }

    std::cout << "[PASS] 3D Taylor-Green vortex decay validated against analytical Navier-Stokes decay rate.\n";
    return 0;
}
