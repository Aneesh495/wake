#include "wake/ispc/lbm_d2q9.isph"
#include "wake/ispc/lbm_d3q19.isph"

// Exported ISPC kernel for 2D D2Q9 Momentum Exchange Method (MEM)
// Integrates net aerodynamic lift, drag, and pitching moment across obstacle boundary links
export void d2q9_compute_momentum_exchange_forces(
    const uniform float f_post_collision[],
    const uniform uint8 cell_type[],
    uniform int nx,
    uniform int ny,
    uniform float center_x,
    uniform float center_y,
    uniform float max_radius,
    uniform float &fx_out,
    uniform float &fy_out,
    uniform float &torque_z_out
) {
    uniform int total_nodes = nx * ny;
    uniform float r_max_sq = max_radius * max_radius;
    varying float sum_fx = 0.0f;
    varying float sum_fy = 0.0f;
    varying float sum_tz = 0.0f;

    for (uniform int q = 1; q < D2Q9_Q; ++q) {
        uniform int cx = d2q9_cx[q];
        uniform int cy = d2q9_cy[q];
        const uniform float* uniform f_q = f_post_collision + q * total_nodes;

        for (uniform int y = 0; y < ny; ++y) {
            uniform int dst_row_offset = y * nx;
            uniform int neighbor_y = y + cy;

            if (neighbor_y >= 0 && neighbor_y < ny) {
                uniform int neighbor_row_offset = neighbor_y * nx;

                foreach (x = 0 ... nx) {
                    varying int idx = dst_row_offset + x;
                    varying uint8 self_type = cell_type[idx];

                    if (self_type == CELL_FLUID) {
                        varying int neighbor_x = (x + cx + nx) % nx;
                        varying int neighbor_idx = neighbor_row_offset + neighbor_x;
                        varying uint8 neighbor_type = cell_type[neighbor_idx];

                        if (neighbor_type == CELL_SOLID || neighbor_type == CELL_MOVING_WALL) {
                            varying float rx = ((float)x + 0.5f * (float)cx) - center_x;
                            varying float ry = ((float)y + 0.5f * (float)cy) - center_y;

                            if (max_radius <= 0.0f || (rx * rx + ry * ry <= r_max_sq)) {
                                varying float f_val = f_q[idx];
                                varying float d_fx = 2.0f * f_val * (float)cx;
                                varying float d_fy = 2.0f * f_val * (float)cy;

                                sum_fx += d_fx;
                                sum_fy += d_fy;
                                sum_tz += (rx * d_fy - ry * d_fx);
                            }
                        }
                    }
                }
            }
        }
    }

    fx_out = reduce_add(sum_fx);
    fy_out = reduce_add(sum_fy);
    torque_z_out = reduce_add(sum_tz);
}

// Exported ISPC kernel for 3D D3Q19 Momentum Exchange Method
export void d3q19_compute_momentum_exchange_forces(
    const uniform float f_post_collision[],
    const uniform uint8 cell_type[],
    uniform int nx,
    uniform int ny,
    uniform int nz,
    uniform float center_x,
    uniform float center_y,
    uniform float center_z,
    uniform float max_radius,
    uniform float &fx_out,
    uniform float &fy_out,
    uniform float &fz_out
) {
    uniform int total_nodes = nx * ny * nz;
    uniform float r_max_sq = max_radius * max_radius;
    varying float sum_fx = 0.0f;
    varying float sum_fy = 0.0f;
    varying float sum_fz = 0.0f;

    for (uniform int q = 1; q < D3Q19_Q; ++q) {
        uniform int cx = d3q19_cx[q];
        uniform int cy = d3q19_cy[q];
        uniform int cz = d3q19_cz[q];
        const uniform float* uniform f_q = f_post_collision + q * total_nodes;

        for (uniform int z = 0; z < nz; ++z) {
            uniform int neighbor_z = z + cz;
            if (neighbor_z < 0 || neighbor_z >= nz) continue;

            for (uniform int y = 0; y < ny; ++y) {
                uniform int neighbor_y = y + cy;
                if (neighbor_y < 0 || neighbor_y >= ny) continue;

                uniform int dst_offset = (z * ny + y) * nx;
                uniform int neighbor_offset = (neighbor_z * ny + neighbor_y) * nx;

                foreach (x = 0 ... nx) {
                    varying int idx = dst_offset + x;
                    varying uint8 self_type = cell_type[idx];

                    if (self_type == CELL_FLUID) {
                        varying int neighbor_x = (x + cx + nx) % nx;
                        varying int neighbor_idx = neighbor_offset + neighbor_x;
                        varying uint8 neighbor_type = cell_type[neighbor_idx];

                        if (neighbor_type == CELL_SOLID || neighbor_type == CELL_MOVING_WALL) {
                            varying float rx = ((float)x + 0.5f * (float)cx) - center_x;
                            varying float ry = ((float)y + 0.5f * (float)cy) - center_y;
                            varying float rz = ((float)z + 0.5f * (float)cz) - center_z;

                            if (max_radius <= 0.0f || (rx * rx + ry * ry + rz * rz <= r_max_sq)) {
                                varying float f_val = f_q[idx];
                                sum_fx += 2.0f * f_val * (float)cx;
                                sum_fy += 2.0f * f_val * (float)cy;
                                sum_fz += 2.0f * f_val * (float)cz;
                            }
                        }
                    }
                }
            }
        }
    }

    fx_out = reduce_add(sum_fx);
    fy_out = reduce_add(sum_fy);
    fz_out = reduce_add(sum_fz);
}
