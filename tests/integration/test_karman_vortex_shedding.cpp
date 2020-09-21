#include "wake/scenarios/karman_vortex.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

using namespace wake;

int main() {
    std::cout << "[INTEGRATION TEST] Running Karman vortex shedding hydrodynamic verification...\n";

    const int NX = 240;
    const int NY = 80;
    const real_t Re = 100.0f;
    const real_t u0 = 0.06f;
    const real_t D = 14.0f;

    KarmanVortexScenario scenario(NX, NY, Re, u0, D);
    scenario.setup();

    real_t initial_mass = scenario.grid().compute_total_mass();
    std::cout << "  Initial domain mass: " << initial_mass << "\n";

    std::vector<real_t> cl_history;
    std::vector<real_t> cd_history;

    const int TOTAL_STEPS = 500;
    for (int step = 1; step <= TOTAL_STEPS; ++step) {
        scenario.step();

        if (std::isnan(scenario.latest_cd())) {
            std::cout << ">>> FIRST NAN AT STEP " << step << " <<<\n";
            for (int y = 0; y < NY; ++y) {
                for (int x = 0; x < NX; ++x) {
                    if (std::isnan(scenario.grid().rho(x, y)) || std::isnan(scenario.grid().ux(x, y))) {
                        std::cout << "NaN detected at cell (" << x << ", " << y << "): type="
                                  << static_cast<int>(scenario.grid().get_cell_type(x, y))
                                  << ", rho=" << scenario.grid().rho(x, y) << ", ux=" << scenario.grid().ux(x, y) << "\n";
                        return 1;
                    }
                }
            }
            return 1;
        }

        cl_history.push_back(scenario.latest_cl());
        cd_history.push_back(scenario.latest_cd());

        if (step % 100 == 0) {
            std::cout << "  Step " << step << "/" << TOTAL_STEPS
                      << " | Cd = " << scenario.latest_cd()
                      << " | Cl = " << scenario.latest_cl()
                      << " | Throughput = " << scenario.current_mlups() << " MLUPS\n";
        }
    }

    std::cout << "Latest forces: Cd=" << scenario.latest_cd() << ", Cl=" << scenario.latest_cl() << "\n";

    real_t final_mass = scenario.grid().compute_total_mass();
    real_t mass_change_pct = std::abs(final_mass - initial_mass) / initial_mass * 100.0f;
    std::cout << "  Final domain mass: " << final_mass << " (Delta: " << mass_change_pct << "%)\n";

    // Physical verification criteria:
    // 1. Drag coefficient must be positive and reasonable
    assert(scenario.latest_cd() > 0.5f);
    assert(scenario.latest_cd() < 50.0f);

    // 2. Maximum fluid velocity must not explode (numerical stability check: u_max < 0.4)
    real_t u_max = scenario.grid().compute_max_speed();
    std::cout << "  Maximum velocity: " << u_max << " (Mach number: " << (u_max / CS) << ")\n";
    assert(u_max < 0.35f);

    // 3. Domain mass must remain conserved within inflow/outflow balance (< 1.5%)
    assert(mass_change_pct < 1.5f);

    std::cout << "[INTEGRATION TEST] Karman vortex shedding verification PASSED!\n";
    return 0;
}
