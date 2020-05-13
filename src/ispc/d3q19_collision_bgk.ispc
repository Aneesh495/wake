#include "wake/ispc/lbm_d3q19.isph"

// Exported ISPC kernel for 3D D3Q19 BGK collision operator with Guo 3D body forcing
export void d3q19_collide_bgk(
    const uniform float f_in[],
    uniform float f_out[],
    const uniform uint8 cell_type[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float uz[],
    uniform float tau,
    uniform float gx,
    uniform float gy,
    uniform float gz,
    uniform int count
) {
    uniform float omega = 1.0f / tau;
    uniform float force_prefactor = (1.0f - 0.5f * omega);

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID) {
            for (uniform int q = 0; q < D3Q19_Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        varying float local_rho = 0.0f;
        varying float jx = 0.0f;
        varying float jy = 0.0f;
        varying float jz = 0.0f;

        varying float f_curr[D3Q19_Q];
        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            varying float fi = f_in[q * count + i];
            f_curr[q] = fi;
            local_rho += fi;
            jx += fi * (float)d3q19_cx[q];
            jy += fi * (float)d3q19_cy[q];
            jz += fi * (float)d3q19_cz[q];
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

        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            varying float ci_x = (float)d3q19_cx[q];
            varying float ci_y = (float)d3q19_cy[q];
            varying float ci_z = (float)d3q19_cz[q];
            varying float ci_u = ci_x * local_ux + ci_y * local_uy + ci_z * local_uz;
            varying float ci_g = ci_x * gx + ci_y * gy + ci_z * gz;

            varying float feq = d3q19_w[q] * local_rho * (
                1.0f + LBM_INV_CS2 * ci_u + LBM_INV_2_CS4 * ci_u * ci_u - LBM_INV_2_CS2 * u_sq
            );

            varying float Fi = force_prefactor * d3q19_w[q] * local_rho * (
                LBM_INV_CS2 * (ci_g - u_dot_g) + 9.0f * ci_u * ci_g
            );

            f_out[q * count + i] = f_curr[q] - omega * (f_curr[q] - feq) + Fi;
        }
    }
}
