#include "wake/ispc/lbm_d3q27.isph"

// Exported ISPC kernel for complete 3D D3Q27 collision with BGK and TRT relaxation
// Eliminates high-order discrete lattice anisotropy artifacts
export void d3q27_collide(
    const uniform float f_in[],
    uniform float f_out[],
    const uniform uint8 cell_type[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float uz[],
    uniform float tau_s,
    uniform float tau_a,
    uniform float gx,
    uniform float gy,
    uniform float gz,
    uniform int count
) {
    uniform float omega_s = 1.0f / tau_s;
    uniform float omega_a = 1.0f / tau_a;
    uniform float prefactor_sym = (1.0f - 0.5f * omega_s);
    uniform float prefactor_asym = (1.0f - 0.5f * omega_a);

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID) {
            for (uniform int q = 0; q < D3Q27_Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        varying float f[D3Q27_Q];
        varying float local_rho = 0.0f;
        varying float jx = 0.0f;
        varying float jy = 0.0f;
        varying float jz = 0.0f;

        for (uniform int q = 0; q < D3Q27_Q; ++q) {
            varying float fi = f_in[q * count + i];
            f[q] = fi;
            local_rho += fi;
            jx += fi * (float)d3q27_cx[q];
            jy += fi * (float)d3q27_cy[q];
            jz += fi * (float)d3q27_cz[q];
        }

        varying float inv_rho = 1.0f / local_rho;
        varying float local_ux = (jx + 0.5f * local_rho * gx) * inv_rho;
        varying float local_uy = (jy + 0.5f * local_rho * gy) * inv_rho;
        varying float local_uz = (jz + 0.5f * local_rho * gz) * inv_rho;

        rho[i] = local_rho;
        ux[i] = local_ux;
        uy[i] = local_uy;
        uz[i] = local_uz;

        varying float u_dot_g = local_ux * gx + local_uy * gy + local_uz * gz;
        varying float u_sq = local_ux * local_ux + local_uy * local_uy + local_uz * local_uz;

        varying float feq[D3Q27_Q];
        for (uniform int q = 0; q < D3Q27_Q; ++q) {
            varying float ci_u = (float)d3q27_cx[q] * local_ux + (float)d3q27_cy[q] * local_uy + (float)d3q27_cz[q] * local_uz;
            feq[q] = d3q27_w[q] * local_rho * (
                1.0f + LBM_INV_CS2 * ci_u + LBM_INV_2_CS4 * ci_u * ci_u - LBM_INV_2_CS2 * u_sq
            );
        }

        for (uniform int q = 0; q < D3Q27_Q; ++q) {
            uniform int inv_q = d3q27_opposite[q];

            varying float f_sym = 0.5f * (f[q] + f[inv_q]);
            varying float f_asym = 0.5f * (f[q] - f[inv_q]);

            varying float feq_sym = 0.5f * (feq[q] + feq[inv_q]);
            varying float feq_asym = 0.5f * (feq[q] - feq[inv_q]);

            varying float ci_x = (float)d3q27_cx[q];
            varying float ci_y = (float)d3q27_cy[q];
            varying float ci_z = (float)d3q27_cz[q];
            varying float ci_u = ci_x * local_ux + ci_y * local_uy + ci_z * local_uz;
            varying float ci_g = ci_x * gx + ci_y * gy + ci_z * gz;

            varying float F_asym = prefactor_asym * d3q27_w[q] * local_rho * (LBM_INV_CS2 * ci_g);
            varying float F_sym  = prefactor_sym  * d3q27_w[q] * local_rho * (-LBM_INV_CS2 * u_dot_g + 9.0f * ci_u * ci_g);
            varying float Fi = F_sym + F_asym;

            f_out[q * count + i] = f[q] - omega_s * (f_sym - feq_sym) - omega_a * (f_asym - feq_asym) + Fi;
        }
    }
}
