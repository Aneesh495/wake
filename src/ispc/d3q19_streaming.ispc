#include "wake/ispc/lbm_d3q19.isph"

// 3D D3Q19 pull-based streaming with periodic boundary wrapping in x, y, and z
export void d3q19_stream_periodic(
    const uniform float f_src[],
    uniform float f_dst[],
    uniform int nx,
    uniform int ny,
    uniform int nz
) {
    uniform int total_nodes = nx * ny * nz;

    for (uniform int q = 0; q < D3Q19_Q; ++q) {
        uniform int cx = d3q19_cx[q];
        uniform int cy = d3q19_cy[q];
        uniform int cz = d3q19_cz[q];

        const uniform float* uniform f_src_q = f_src + q * total_nodes;
        uniform float* uniform f_dst_q = f_dst + q * total_nodes;

        for (uniform int z = 0; z < nz; ++z) {
            uniform int src_z = (z - cz + nz) % nz;
            for (uniform int y = 0; y < ny; ++y) {
                uniform int src_y = (y - cy + ny) % ny;
                uniform int dst_offset = (z * ny + y) * nx;
                uniform int src_offset = (src_z * ny + src_y) * nx;

                foreach (x = 0 ... nx) {
                    varying int src_x = (x - cx + nx) % nx;
                    f_dst_q[dst_offset + x] = f_src_q[src_offset + src_x];
                }
            }
        }
    }
}

// 3D D3Q19 streaming with integrated halfway bounce-back on solid obstacles
// Supports x-periodic or channel flow with solid walls in y and z
export void d3q19_stream_with_boundaries(
    const uniform float f_src[],
    uniform float f_dst[],
    const uniform uint8 cell_type[],
    uniform int nx,
    uniform int ny,
    uniform int nz,
    uniform float uw_x,
    uniform float uw_y,
    uniform float uw_z
) {
    uniform int total_nodes = nx * ny * nz;

    varying float momentum_term[D3Q19_Q];
    for (uniform int q = 0; q < D3Q19_Q; ++q) {
        varying float ci_uw = (float)d3q19_cx[q] * uw_x + (float)d3q19_cy[q] * uw_y + (float)d3q19_cz[q] * uw_z;
        momentum_term[q] = 2.0f * d3q19_w[q] * LBM_RHO_0 * LBM_INV_CS2 * ci_uw;
    }

    for (uniform int q = 0; q < D3Q19_Q; ++q) {
        uniform int cx = d3q19_cx[q];
        uniform int cy = d3q19_cy[q];
        uniform int cz = d3q19_cz[q];
        uniform int opp_q = d3q19_opposite[q];

        const uniform float* uniform f_src_q = f_src + q * total_nodes;
        const uniform float* uniform f_src_opp = f_src + opp_q * total_nodes;
        uniform float* uniform f_dst_q = f_dst + q * total_nodes;

        for (uniform int z = 0; z < nz; ++z) {
            uniform int src_z = z - cz;
            for (uniform int y = 0; y < ny; ++y) {
                uniform int src_y = y - cy;
                uniform int dst_offset = (z * ny + y) * nx;

                foreach (x = 0 ... nx) {
                    varying int idx = dst_offset + x;
                    varying uint8 self_type = cell_type[idx];

                    if (self_type == CELL_SOLID) {
                        f_dst_q[idx] = f_src_q[idx];
                        continue;
                    }

                    varying int src_x = (x - cx + nx) % nx; // periodic in x

                    if (src_y >= 0 && src_y < ny && src_z >= 0 && src_z < nz) {
                        varying int src_idx = (src_z * ny + src_y) * nx + src_x;
                        varying uint8 neighbor_type = cell_type[src_idx];

                        if (neighbor_type == CELL_SOLID) {
                            f_dst_q[idx] = f_src_opp[idx];
                        } else if (neighbor_type == CELL_MOVING_WALL) {
                            f_dst_q[idx] = f_src_opp[idx] + momentum_term[q];
                        } else {
                            f_dst_q[idx] = f_src_q[src_idx];
                        }
                    } else {
                        // Halfway bounce-back at domain domain wall boundaries
                        f_dst_q[idx] = f_src_opp[idx];
                    }
                }
            }
        }
    }
}
