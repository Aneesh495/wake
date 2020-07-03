#include "wake/ispc/lbm_d2q9.isph"

// Exported ISPC kernel for thermal D2Q9 collision and Boussinesq buoyancy coupling
export void thermal_collide_d2q9(
    const uniform float g_in[],
    uniform float g_out[],
    const uniform uint8 cell_type[],
    const uniform float ux[],
    const uniform float uy[],
    uniform float temperature[],
    uniform float tau_T,
    uniform int count
) {
    uniform float omega_T = 1.0f / tau_T;

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID) {
            for (uniform int q = 0; q < D2Q9_Q; ++q) {
                g_out[q * count + i] = g_in[q * count + i];
            }
            continue;
        }

        varying float local_T = 0.0f;
        varying float g_curr[D2Q9_Q];

        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            varying float gi = g_in[q * count + i];
            g_curr[q] = gi;
            local_T += gi;
        }

        temperature[i] = local_T;
        varying float vx = ux[i];
        varying float vy = uy[i];

        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            varying float ci_u = (float)d2q9_cx[q] * vx + (float)d2q9_cy[q] * vy;
            // First-order thermal equilibrium: geq = w_i * T * (1 + 3 * c_i . u)
            varying float geq = d2q9_w[q] * local_T * (1.0f + LBM_INV_CS2 * ci_u);
            g_out[q * count + i] = g_curr[q] - omega_T * (g_curr[q] - geq);
        }
    }
}

// Boussinesq thermal buoyancy body force computation
// Adds buoyancy acceleration: gy_buoy = beta * (T - T_ref) * g_gravity
export void thermal_compute_buoyancy_force_2d(
    const uniform float temperature[],
    uniform float gy_body[],
    uniform float t_ref,
    uniform float beta,
    uniform float g_gravity,
    uniform int count
) {
    foreach (i = 0 ... count) {
        varying float T = temperature[i];
        gy_body[i] += beta * (T - t_ref) * g_gravity;
    }
}

// Dirichlet isothermal wall boundary condition: g_opp = -g_in + 2 * w_i * T_wall
export void thermal_apply_isothermal_wall(
    uniform float g[],
    uniform int nx,
    uniform int ny,
    uniform int y_wall,
    uniform float t_wall
) {
    uniform int total_nodes = nx * ny;

    foreach (x = 0 ... nx) {
        varying int idx = y_wall * nx + x;

        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            uniform int opp_q = d2q9_opposite[q];
            varying float g_val = g[q * total_nodes + idx];
            g[opp_q * total_nodes + idx] = -g_val + 2.0f * d2q9_w[q] * t_wall;
        }
    }
}
