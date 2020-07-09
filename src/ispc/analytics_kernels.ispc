#include "wake/ispc/lbm_common.isph"

// Exported ISPC kernel for 2D vorticity computation: omega_z = duy/dx - dux/dy
export void compute_vorticity_2d(
    const uniform float ux[],
    const uniform float uy[],
    uniform float vorticity[],
    uniform int nx,
    uniform int ny
) {
    for (uniform int y = 0; y < ny; ++y) {
        uniform int ym = (y - 1 + ny) % ny;
        uniform int yp = (y + 1) % ny;
        uniform int row_curr = y * nx;
        uniform int row_prev = ym * nx;
        uniform int row_next = yp * nx;

        foreach (x = 0 ... nx) {
            varying int xm = (x - 1 + nx) % nx;
            varying int xp = (x + 1) % nx;

            // Central differences
            varying float duy_dx = 0.5f * (uy[row_curr + xp] - uy[row_curr + xm]);
            varying float dux_dy = 0.5f * (ux[row_next + x] - ux[row_prev + x]);

            vorticity[row_curr + x] = duy_dx - dux_dy;
        }
    }
}

// Exported ISPC kernel for 3D Vorticity vector and Q-criterion computation
// Q = 0.5 * (||Omega||^2 - ||S||^2)
export void compute_vorticity_and_q_criterion_3d(
    const uniform float ux[],
    const uniform float uy[],
    const uniform float uz[],
    uniform float omega_x[],
    uniform float omega_y[],
    uniform float omega_z[],
    uniform float q_criterion[],
    uniform int nx,
    uniform int ny,
    uniform int nz
) {
    for (uniform int z = 0; z < nz; ++z) {
        uniform int zm = (z - 1 + nz) % nz;
        uniform int zp = (z + 1) % nz;

        for (uniform int y = 0; y < ny; ++y) {
            uniform int ym = (y - 1 + ny) % ny;
            uniform int yp = (y + 1) % ny;

            uniform int idx_center = (z * ny + y) * nx;
            uniform int idx_yp = (z * ny + yp) * nx;
            uniform int idx_ym = (z * ny + ym) * nx;
            uniform int idx_zp = (zp * ny + y) * nx;
            uniform int idx_zm = (zm * ny + y) * nx;

            foreach (x = 0 ... nx) {
                varying int xm = (x - 1 + nx) % nx;
                varying int xp = (x + 1) % nx;

                varying int curr = idx_center + x;

                // 9 components of velocity gradient tensor grad(u)
                varying float dux_dx = 0.5f * (ux[idx_center + xp] - ux[idx_center + xm]);
                varying float dux_dy = 0.5f * (ux[idx_yp + x] - ux[idx_ym + x]);
                varying float dux_dz = 0.5f * (ux[idx_zp + x] - ux[idx_zm + x]);

                varying float duy_dx = 0.5f * (uy[idx_center + xp] - uy[idx_center + xm]);
                varying float duy_dy = 0.5f * (uy[idx_yp + x] - uy[idx_ym + x]);
                varying float duy_dz = 0.5f * (uy[idx_zp + x] - uy[idx_zm + x]);

                varying float duz_dx = 0.5f * (uz[idx_center + xp] - uz[idx_center + xm]);
                varying float duz_dy = 0.5f * (uz[idx_yp + x] - uz[idx_ym + x]);
                varying float duz_dz = 0.5f * (uz[idx_zp + x] - uz[idx_zm + x]);

                // Vorticity components: curl(u)
                varying float wx = duz_dy - duy_dz;
                varying float wy = dux_dz - duz_dx;
                varying float wz = duy_dx - dux_dy;

                omega_x[curr] = wx;
                omega_y[curr] = wy;
                omega_z[curr] = wz;

                // Symmetric strain rate tensor S = 0.5 * (grad(u) + grad(u)^T)
                varying float s_xx = dux_dx;
                varying float s_yy = duy_dy;
                varying float s_zz = duz_dz;
                varying float s_xy = 0.5f * (dux_dy + duy_dx);
                varying float s_yz = 0.5f * (duy_dz + duz_dy);
                varying float s_xz = 0.5f * (dux_dz + duz_dx);

                varying float norm_s_sq = (s_xx * s_xx + s_yy * s_yy + s_zz * s_zz)
                                        + 2.0f * (s_xy * s_xy + s_yz * s_yz + s_xz * s_xz);

                // Antisymmetric rotation tensor Omega = 0.5 * (grad(u) - grad(u)^T)
                varying float o_xy = 0.5f * (dux_dy - duy_dx);
                varying float o_yz = 0.5f * (duy_dz - duz_dy);
                varying float o_xz = 0.5f * (dux_dz - duz_dx);

                varying float norm_o_sq = 2.0f * (o_xy * o_xy + o_yz * o_yz + o_xz * o_xz);

                // Q-criterion: Q = 0.5 * (||Omega||^2 - ||S||^2)
                q_criterion[curr] = 0.5f * (norm_o_sq - norm_s_sq);
            }
        }
    }
}
