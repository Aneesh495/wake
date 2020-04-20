#include "wake/ispc/lbm_d2q9.isph"

// Pull-based streaming kernel for D2Q9 lattice with periodic wrapping
// Row-wise uniform y decomposition with SIMD vectorization across x columns
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

// Pull-based streaming for general domains with interior propagation
export void d2q9_stream_interior(
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
            uniform int src_y = y - cy;
            uniform int dst_row_offset = y * nx;

            if (src_y >= 0 && src_y < ny) {
                uniform int src_row_offset = src_y * nx;
                foreach (x = 0 ... nx) {
                    varying int src_x = x - cx;
                    if (src_x >= 0 && src_x < nx) {
                        f_dst_q[dst_row_offset + x] = f_src_q[src_row_offset + src_x];
                    } else {
                        f_dst_q[dst_row_offset + x] = f_src_q[dst_row_offset + x];
                    }
                }
            } else {
                foreach (x = 0 ... nx) {
                    f_dst_q[dst_row_offset + x] = f_src_q[dst_row_offset + x];
                }
            }
        }
    }
}
