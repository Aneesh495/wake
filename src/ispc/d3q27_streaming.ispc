#include "wake/ispc/lbm_d3q27.isph"

// 3D D3Q27 pull-based streaming with 3D periodic boundary wrapping
export void d3q27_stream_periodic(
    const uniform float f_src[],
    uniform float f_dst[],
    uniform int nx,
    uniform int ny,
    uniform int nz
) {
    uniform int total_nodes = nx * ny * nz;

    for (uniform int q = 0; q < D3Q27_Q; ++q) {
        uniform int cx = d3q27_cx[q];
        uniform int cy = d3q27_cy[q];
        uniform int cz = d3q27_cz[q];

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
