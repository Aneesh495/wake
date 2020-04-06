#include "wake/ispc/lbm_d2q9.isph"

// Exported ISPC kernel for D2Q9 Bhatnagar-Gross-Krook (BGK) collision
// Supports Guo body-forcing for external gravity and pressure gradients
export void d2q9_collide_bgk(
    const uniform float f_in[],
    uniform float f_out[],
    const uniform uint8 cell_type[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float tau,
    uniform float gx,
    uniform float gy,
    uniform int count
) {
    uniform float omega = 1.0f / tau;
    uniform float force_prefactor = (1.0f - 0.5f * omega);

    foreach (i = 0 ... count) {
        // Skip non-fluid nodes (boundary nodes handled separately)
        if (cell_type[i] == CELL_SOLID) {
            // Keep distributions intact for solid nodes
            for (uniform int q = 0; q < D2Q9_Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        // Gather 9 distributions
        varying float f0 = f_in[0 * count + i];
        varying float f1 = f_in[1 * count + i];
        varying float f2 = f_in[2 * count + i];
        varying float f3 = f_in[3 * count + i];
        varying float f4 = f_in[4 * count + i];
        varying float f5 = f_in[5 * count + i];
        varying float f6 = f_in[6 * count + i];
        varying float f7 = f_in[7 * count + i];
        varying float f8 = f_in[8 * count + i];

        // Density
        varying float local_rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
        varying float inv_rho = 1.0f / local_rho;

        // Momentum
        varying float jx = (f1 - f3) + (f5 - f6 - f7 + f8);
        varying float jy = (f2 - f4) + (f5 + f6 - f7 - f8);

        // Guo velocity shift due to body force
        varying float local_ux = (jx + 0.5f * local_rho * gx) * inv_rho;
        varying float local_uy = (jy + 0.5f * local_rho * gy) * inv_rho;

        rho[i] = local_rho;
        ux[i] = local_ux;
        uy[i] = local_uy;

        varying float u_dot_g = local_ux * gx + local_uy * gy;
        varying float u_sq = local_ux * local_ux + local_uy * local_uy;

        // Apply BGK relaxation with Guo forcing term across all 9 directions
        varying float f_curr[D2Q9_Q];
        f_curr[0] = f0; f_curr[1] = f1; f_curr[2] = f2;
        f_curr[3] = f3; f_curr[4] = f4; f_curr[5] = f5;
        f_curr[6] = f6; f_curr[7] = f7; f_curr[8] = f8;

        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            varying float ci_x = (float)d2q9_cx[q];
            varying float ci_y = (float)d2q9_cy[q];
            varying float ci_u = ci_x * local_ux + ci_y * local_uy;
            varying float ci_g = ci_x * gx + ci_y * gy;

            // Equilibrium distribution
            varying float feq = d2q9_w[q] * local_rho * (
                1.0f + LBM_INV_CS2 * ci_u + LBM_INV_2_CS4 * ci_u * ci_u - LBM_INV_2_CS2 * u_sq
            );

            // Guo body force source term
            varying float Fi = force_prefactor * d2q9_w[q] * local_rho * (
                LBM_INV_CS2 * (ci_g - u_dot_g) + 9.0f * ci_u * ci_g
            );

            // Post-collision distribution
            f_out[q * count + i] = f_curr[q] - omega * (f_curr[q] - feq) + Fi;
        }
    }
}
