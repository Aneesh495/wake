#include "wake/ispc/lbm_d2q9.isph"

// Exported ISPC kernel for D2Q9 Two-Relaxation-Time (TRT) collision operator
// Decouples symmetric (viscous) and anti-symmetric (advective) relaxation rates
// Magic parameter Lambda = (tau_s - 0.5) * (tau_a - 0.5) stabilizes boundary locations
// Exact TRT Guo body force decomposition: symmetric (omega_s) and anti-symmetric (omega_a)
export void d2q9_collide_trt(
    const uniform float f_in[],
    uniform float f_out[],
    const uniform uint8 cell_type[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float tau_s,
    uniform float tau_a,
    uniform float gx,
    uniform float gy,
    uniform int count
) {
    uniform float omega_s = 1.0f / tau_s;
    uniform float omega_a = 1.0f / tau_a;
    uniform float prefactor_sym = (1.0f - 0.5f * omega_s);
    uniform float prefactor_asym = (1.0f - 0.5f * omega_a);

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID) {
            for (uniform int q = 0; q < D2Q9_Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        varying float f[D2Q9_Q];
        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            f[q] = f_in[q * count + i];
        }

        // Density
        varying float local_rho = 0.0f;
        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            local_rho += f[q];
        }
        varying float inv_rho = 1.0f / local_rho;

        // Momentum
        varying float jx = (f[1] - f[3]) + (f[5] - f[6] - f[7] + f[8]);
        varying float jy = (f[2] - f[4]) + (f[5] + f[6] - f[7] - f[8]);

        // Guo velocity shift
        varying float local_ux = (jx + 0.5f * local_rho * gx) * inv_rho;
        varying float local_uy = (jy + 0.5f * local_rho * gy) * inv_rho;

        rho[i] = local_rho;
        ux[i] = local_ux;
        uy[i] = local_uy;

        varying float u_dot_g = local_ux * gx + local_uy * gy;
        varying float u_sq = local_ux * local_ux + local_uy * local_uy;

        // Compute equilibrium distributions
        varying float feq[D2Q9_Q];
        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            varying float ci_u = (float)d2q9_cx[q] * local_ux + (float)d2q9_cy[q] * local_uy;
            feq[q] = d2q9_w[q] * local_rho * (
                1.0f + LBM_INV_CS2 * ci_u + LBM_INV_2_CS4 * ci_u * ci_u - LBM_INV_2_CS2 * u_sq
            );
        }

        // TRT relaxation with symmetric and anti-symmetric decomposition
        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            uniform int inv_q = d2q9_opposite[q];

            varying float f_sym = 0.5f * (f[q] + f[inv_q]);
            varying float f_asym = 0.5f * (f[q] - f[inv_q]);

            varying float feq_sym = 0.5f * (feq[q] + feq[inv_q]);
            varying float feq_asym = 0.5f * (feq[q] - feq[inv_q]);

            varying float ci_x = (float)d2q9_cx[q];
            varying float ci_y = (float)d2q9_cy[q];
            varying float ci_u = ci_x * local_ux + ci_y * local_uy;
            varying float ci_g = ci_x * gx + ci_y * gy;

            // Decomposed Guo force: odd part (ci_g) relaxes with omega_a; even part relaxes with omega_s
            varying float F_asym = prefactor_asym * d2q9_w[q] * local_rho * (LBM_INV_CS2 * ci_g);
            varying float F_sym  = prefactor_sym  * d2q9_w[q] * local_rho * (-LBM_INV_CS2 * u_dot_g + 9.0f * ci_u * ci_g);
            varying float Fi = F_sym + F_asym;

            f_out[q * count + i] = f[q] - omega_s * (f_sym - feq_sym) - omega_a * (f_asym - feq_asym) + Fi;
        }
    }
}
