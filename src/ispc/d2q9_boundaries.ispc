#include "wake/ispc/lbm_d2q9.isph"

// Half-way bounce-back kernel for solid obstacle nodes
// If wall moves with velocity (uw_x, uw_y), momentum transfer term is added
export void d2q9_apply_bounce_back(
    uniform float f[],
    const uniform uint8 cell_type[],
    uniform float uw_x,
    uniform float uw_y,
    uniform int count
) {
    varying float momentum_term[D2Q9_Q];
    for (uniform int q = 0; q < D2Q9_Q; ++q) {
        varying float ci_uw = (float)d2q9_cx[q] * uw_x + (float)d2q9_cy[q] * uw_y;
        momentum_term[q] = 2.0f * d2q9_w[q] * LBM_RHO_0 * LBM_INV_CS2 * ci_uw;
    }

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID || cell_type[i] == CELL_MOVING_WALL) {
            // Read all current distributions
            varying float f_curr[D2Q9_Q];
            for (uniform int q = 0; q < D2Q9_Q; ++q) {
                f_curr[q] = f[q * count + i];
            }

            // Swap opposite directions
            // f[inv_q] = f_curr[q] - momentum_term[q]
            f[0 * count + i] = f_curr[0];
            f[1 * count + i] = f_curr[3] + momentum_term[1];
            f[2 * count + i] = f_curr[4] + momentum_term[2];
            f[3 * count + i] = f_curr[1] + momentum_term[3];
            f[4 * count + i] = f_curr[2] + momentum_term[4];
            f[5 * count + i] = f_curr[7] + momentum_term[5];
            f[6 * count + i] = f_curr[8] + momentum_term[6];
            f[7 * count + i] = f_curr[5] + momentum_term[7];
            f[8 * count + i] = f_curr[6] + momentum_term[8];
        }
    }
}

// Zou-He velocity inflow boundary condition at west wall (x = 0)
// Prescribes velocity (u_in_x, u_in_y) and computes incoming unknown distributions f1, f5, f8
export void d2q9_apply_zou_he_inflow_west(
    uniform float f[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float u_in_x,
    uniform float u_in_y,
    uniform int nx,
    uniform int ny
) {
    uniform int total_nodes = nx * ny;

    foreach (y = 1 ... ny - 1) {
        varying int idx = y * nx + 0;

        varying float f0 = f[0 * total_nodes + idx];
        varying float f2 = f[2 * total_nodes + idx];
        varying float f3 = f[3 * total_nodes + idx];
        varying float f4 = f[4 * total_nodes + idx];
        varying float f6 = f[6 * total_nodes + idx];
        varying float f7 = f[7 * total_nodes + idx];

        // Density derived from velocity and known distributions
        varying float rho_in = (f0 + f2 + f4 + 2.0f * (f3 + f6 + f7)) / (1.0f - u_in_x);

        // Unknown incoming populations
        varying float f1 = f3 + (2.0f / 3.0f) * rho_in * u_in_x;
        varying float f5 = f7 - 0.5f * (f2 - f4) + (1.0f / 6.0f) * rho_in * u_in_x + 0.5f * rho_in * u_in_y;
        varying float f8 = f6 + 0.5f * (f2 - f4) + (1.0f / 6.0f) * rho_in * u_in_x - 0.5f * rho_in * u_in_y;

        f[1 * total_nodes + idx] = f1;
        f[5 * total_nodes + idx] = f5;
        f[8 * total_nodes + idx] = f8;

        rho[idx] = rho_in;
        ux[idx] = u_in_x;
        uy[idx] = u_in_y;
    }
}

// Zou-He pressure outflow boundary condition at east wall (x = nx - 1)
// Prescribes density rho_out and computes incoming unknown distributions f3, f6, f7
export void d2q9_apply_zou_he_outflow_east(
    uniform float f[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float rho_out,
    uniform int nx,
    uniform int ny
) {
    uniform int total_nodes = nx * ny;

    foreach (y = 1 ... ny - 1) {
        varying int idx = y * nx + (nx - 1);

        varying float f0 = f[0 * total_nodes + idx];
        varying float f1 = f[1 * total_nodes + idx];
        varying float f2 = f[2 * total_nodes + idx];
        varying float f4 = f[4 * total_nodes + idx];
        varying float f5 = f[5 * total_nodes + idx];
        varying float f8 = f[8 * total_nodes + idx];

        varying float u_out_x = -1.0f + (f0 + f2 + f4 + 2.0f * (f1 + f5 + f8)) / rho_out;
        varying float u_out_y = 0.0f;

        varying float f3 = f1 - (2.0f / 3.0f) * rho_out * u_out_x;
        varying float f6 = f8 - 0.5f * (f2 - f4) - (1.0f / 6.0f) * rho_out * u_out_x + 0.5f * rho_out * u_out_y;
        varying float f7 = f5 + 0.5f * (f2 - f4) - (1.0f / 6.0f) * rho_out * u_out_x - 0.5f * rho_out * u_out_y;

        f[3 * total_nodes + idx] = f3;
        f[6 * total_nodes + idx] = f6;
        f[7 * total_nodes + idx] = f7;

        rho[idx] = rho_out;
        ux[idx] = u_out_x;
        uy[idx] = u_out_y;
    }
}

// Convective non-reflecting outflow boundary condition at east wall (x = nx - 1)
// Discretizes df/dt + U_conv * df/dx = 0
export void d2q9_apply_convective_outflow_east(
    uniform float f[],
    const uniform float f_prev[],
    uniform float u_conv,
    uniform int nx,
    uniform int ny
) {
    uniform int total_nodes = nx * ny;
    uniform float alpha = u_conv / (1.0f + u_conv);
    uniform float beta = 1.0f / (1.0f + u_conv);

    foreach (y = 1 ... ny - 1) {
        varying int idx_out = y * nx + (nx - 1);
        varying int idx_in  = y * nx + (nx - 2);

        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            varying float f_prev_out = f_prev[q * total_nodes + idx_out];
            varying float f_curr_in  = f[q * total_nodes + idx_in];
            f[q * total_nodes + idx_out] = beta * f_prev_out + alpha * f_curr_in;
        }
    }
}
