#include "wake/grid/lattice_grid_2d.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace wake;

int main() {
    std::cout << "[VALIDATION] Running 2D Poiseuille channel flow analytical convergence test...\n";

    const int NX = 8;
    const int NY = 33; // y=0 and y=32 are solid walls. Fluid nodes y=1..31.
    const real_t tau = 0.8f;
    const real_t nu = (1.0f / 3.0f) * (tau - 0.5f); // 0.1
    const real_t gx = 1e-4f; // Body force
    const real_t gy = 0.0f;
    const real_t H_eff = static_cast<real_t>(NY - 2); // Effective width between half-way walls: 31.0
    const real_t u_max_analytical = (gx * H_eff * H_eff) / (8.0f * nu);

    LatticeGrid2D grid(NX, NY, tau);

    // Set top and bottom no-slip solid walls
    for (int x = 0; x < NX; ++x) {
        grid.set_cell_type(x, 0, CellType::Solid);
        grid.set_cell_type(x, NY - 1, CellType::Solid);
    }

    // Run simulation to full steady state (tau_diff ~ 9600 steps)
    const int max_steps = 7000;
    for (int step = 0; step < max_steps; ++step) {
        grid.step(CollisionType::TRT, gx, gy);
    }

    // Extract velocity profile across fluid nodes y=1..NY-2
    real_t l2_error_numerator = 0.0f;
    real_t l2_error_denominator = 0.0f;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  y | Simulated u_x | Analytical u_x | Error\n";
    std::cout << "----+---------------+----------------+----------\n";

    for (int y = 1; y < NY - 1; ++y) {
        real_t u_sim = 0.0f;
        for (int x = 0; x < NX; ++x) {
            u_sim += grid.ux(x, y);
        }
        u_sim /= static_cast<real_t>(NX);

        // Distance from bottom halfway wall at y = 0.5
        real_t dist_from_wall = static_cast<real_t>(y) - 0.5f;
        real_t u_ana = (gx / (2.0f * nu)) * dist_from_wall * (H_eff - dist_from_wall);

        real_t diff = std::abs(u_sim - u_ana);
        l2_error_numerator += diff * diff;
        l2_error_denominator += u_ana * u_ana;

        if (y % 4 == 1) {
            std::cout << std::setw(3) << y << " | "
                      << std::setw(13) << u_sim << " | "
                      << std::setw(14) << u_ana << " | "
                      << std::setw(8) << diff << "\n";
        }
    }

    real_t relative_l2_error = std::sqrt(l2_error_numerator / l2_error_denominator);
    std::cout << "Relative L2 Error: " << relative_l2_error * 100.0f << " %\n";
    std::cout << "Analytical u_max: " << u_max_analytical << ", Simulated u_max: " << grid.compute_max_speed() << "\n";

    if (relative_l2_error > 0.015f) { // Must be within 1.5%
        std::cerr << "Poiseuille validation failed! L2 error too high: " << relative_l2_error << "\n";
        return 1;
    }

    std::cout << "[PASS] Poiseuille channel flow validation verified against analytical solution.\n";
    return 0;
}
