#include "wake/ispc/lbm_d2q9.isph"

// Pull-based streaming kernel for D2Q9 lattice with periodic wrapping
export void d2q9_stream_periodic(
    const uniform float f_src[],
    uniform float f_dst[],
    uniform int nx,
    uniform int ny
) {
    uniform int total_nodes = nx * ny;

    for (uniform int q = 0; q < D2Q9_Q; ++q) {
        uniform int cx = d2q9_cx[q];
        uniform int cy = d2q9_cy[q];
        const uniform float* uniform f_src_q = f_src + q * total_nodes;
        uniform float* uniform f_dst_q = f_dst + q * total_nodes;

        for (uniform int y = 0; y < ny; ++y) {
            uniform int src_y = (y - cy + ny) % ny;
            uniform int src_row_offset = src_y * nx;
            uniform int dst_row_offset = y * nx;

            foreach (x = 0 ... nx) {
                varying int src_x = (x - cx + nx) % nx;
                f_dst_q[dst_row_offset + x] = f_src_q[src_row_offset + src_x];
            }
        }
    }
}

// Pull-based streaming with integrated half-way bounce-back on solid boundary links
// For any link pointing into a solid node, reflects the opposite distribution from the current fluid node
export void d2q9_stream_with_boundaries(
    const uniform float f_src[],
    uniform float f_dst[],
    const uniform uint8 cell_type[],
    uniform int nx,
    uniform int ny,
    uniform float uw_x,
    uniform float uw_y
) {
    uniform int total_nodes = nx * ny;

    varying float momentum_term[D2Q9_Q];
    for (uniform int q = 0; q < D2Q9_Q; ++q) {
        varying float ci_uw = (float)d2q9_cx[q] * uw_x + (float)d2q9_cy[q] * uw_y;
        momentum_term[q] = 2.0f * d2q9_w[q] * LBM_RHO_0 * LBM_INV_CS2 * ci_uw;
    }

    for (uniform int q = 0; q < D2Q9_Q; ++q) {
        uniform int cx = d2q9_cx[q];
        uniform int cy = d2q9_cy[q];
        uniform int opp_q = d2q9_opposite[q];

        const uniform float* uniform f_src_q = f_src + q * total_nodes;
        const uniform float* uniform f_src_opp = f_src + opp_q * total_nodes;
        uniform float* uniform f_dst_q = f_dst + q * total_nodes;

        for (uniform int y = 0; y < ny; ++y) {
            uniform int dst_row_offset = y * nx;
            uniform int src_y = y - cy;

            foreach (x = 0 ... nx) {
                varying int idx = dst_row_offset + x;
                varying uint8 self_type = cell_type[idx];

                if (self_type == CELL_SOLID) {
                    // Solid nodes do not participate in fluid transport
                    f_dst_q[idx] = f_src_q[idx];
                    continue;
                }

                varying int src_x = (x - cx + nx) % nx; // periodic along x

                if (src_y >= 0 && src_y < ny) {
                    varying int src_idx = src_y * nx + src_x;
                    varying uint8 neighbor_type = cell_type[src_idx];

                    if (neighbor_type == CELL_SOLID) {
                        // Halfway bounce-back: incoming q receives current node's opposite q distribution
                        f_dst_q[idx] = f_src_opp[idx];
                    } else if (neighbor_type == CELL_MOVING_WALL) {
                        f_dst_q[idx] = f_src_opp[idx] + momentum_term[q];
                    } else {
                        // Normal fluid streaming
                        f_dst_q[idx] = f_src_q[src_idx];
                    }
                } else {
                    // Domain top or bottom boundary without wall node (treat as no-slip halfway wall)
                    f_dst_q[idx] = f_src_opp[idx];
                }
            }
        }
    }
}
