#include "wake/ispc/lbm_d3q19.isph"

// Exported ISPC kernel for Bouzidi-Firdaouss-Lallemand (BFL) interpolated bounce-back
// on curved 3D obstacle boundaries.
// Replaces first-order stair-stepped approximation with second-order spatial interpolation
// based on intersection fraction delta = |x_f - x_w| / |x_f - x_s|
export void d3q19_apply_curved_bfl(
    uniform float f_dst[],
    const uniform float f_src[],
    const uniform int fluid_indices[],
    const uniform int q_dirs[],
    const uniform float delta_fractions[],
    const uniform float uw_x[],
    const uniform float uw_y[],
    const uniform float uw_z[],
    uniform int nx,
    uniform int ny,
    uniform int nz,
    uniform int num_boundary_links
) {
    uniform int total_nodes = nx * ny * nz;

    foreach (link = 0 ... num_boundary_links) {
        varying int idx_f = fluid_indices[link];
        varying int q = q_dirs[link];
        varying float delta = delta_fractions[link];

        varying int opp_q = d3q19_opposite[q];
        varying int cx = d3q19_cx[q];
        varying int cy = d3q19_cy[q];
        varying int cz = d3q19_cz[q];

        // Wall velocity
        varying float vx_w = uw_x[link];
        varying float vy_w = uw_y[link];
        varying float vz_w = uw_z[link];
        varying float ci_uw = (float)cx * vx_w + (float)cy * vy_w + (float)cz * vz_w;

        varying float f_i_xf = f_src[q * total_nodes + idx_f];
        varying float f_opp_xf = f_src[opp_q * total_nodes + idx_f];

        varying float reflected_val;

        if (delta < 0.5f) {
            // Fluid neighbor further upstream: x_ff = x_f - c_i
            varying int z_f = idx_f / (nx * ny);
            varying int rem = idx_f - z_f * (nx * ny);
            varying int y_f = rem / nx;
            varying int x_f = rem - y_f * nx;

            varying int x_ff = (x_f - cx + nx) % nx;
            varying int y_ff = y_f - cy;
            varying int z_ff = z_f - cz;

            varying float f_i_xff = f_i_xf;
            if (y_ff >= 0 && y_ff < ny && z_ff >= 0 && z_ff < nz) {
                varying int idx_ff = (z_ff * ny + y_ff) * nx + x_ff;
                f_i_xff = f_src[q * total_nodes + idx_ff];
            }

            varying float mom = 2.0f * d3q19_w[q] * LBM_RHO_0 * LBM_INV_CS2 * ci_uw;
            reflected_val = 2.0f * delta * f_i_xf + (1.0f - 2.0f * delta) * f_i_xff - mom;
        } else {
            varying float inv_2delta = 1.0f / (2.0f * delta);
            varying float mom = (1.0f / delta) * d3q19_w[q] * LBM_RHO_0 * LBM_INV_CS2 * ci_uw;
            reflected_val = inv_2delta * f_i_xf + (1.0f - inv_2delta) * f_opp_xf - mom;
        }

        f_dst[opp_q * total_nodes + idx_f] = reflected_val;
    }
}
