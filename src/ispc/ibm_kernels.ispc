#include "wake/ispc/lbm_common.isph"

// Roma-Peskin regularized 4-point discrete delta function component
static inline varying float peskin_kernel_1d(varying float r) {
    varying float abs_r = abs(r);
    if (abs_r >= 2.0f) return 0.0f;
    return 0.25f * (1.0f + cos(0.5f * 3.14159265358979323846f * abs_r));
}

// 2D Immersed Boundary velocity interpolation
// Interpolates fluid velocity from Eulerian grid to Lagrangian markers
export void ibm_interpolate_velocity_2d(
    const uniform float ux[],
    const uniform float uy[],
    const uniform float marker_x[],
    const uniform float marker_y[],
    uniform float marker_ux[],
    uniform float marker_uy[],
    uniform int nx,
    uniform int ny,
    uniform int num_markers
) {
    foreach (m = 0 ... num_markers) {
        varying float mx = marker_x[m];
        varying float my = marker_y[m];

        varying float interp_ux = 0.0f;
        varying float interp_uy = 0.0f;

        for (uniform int iy = -2; iy <= 2; ++iy) {
            varying int py = (int)floor(my) + iy;
            if (py >= 0 && py < ny) {
                varying float wy = peskin_kernel_1d((float)py - my);
                varying int row_offset = py * nx;

                for (uniform int ix = -2; ix <= 2; ++ix) {
                    varying int px = (int)floor(mx) + ix;
                    varying int px_wrapped = (px % nx + nx) % nx; // periodic along x
                    varying float wx = peskin_kernel_1d((float)px - mx);
                    varying float weight = wx * wy;

                    varying int euler_idx = row_offset + px_wrapped;
                    interp_ux += weight * ux[euler_idx];
                    interp_uy += weight * uy[euler_idx];
                }
            }
        }

        marker_ux[m] = interp_ux;
        marker_uy[m] = interp_uy;
    }
}

// 2D Immersed Boundary force spreading
// Distributes Lagrangian restoring forces to Eulerian lattice body force array (gx, gy)
export void ibm_spread_force_2d(
    const uniform float marker_x[],
    const uniform float marker_y[],
    const uniform float marker_fx[],
    const uniform float marker_fy[],
    const uniform float marker_ds[],
    uniform float euler_gx[],
    uniform float euler_gy[],
    uniform int nx,
    uniform int ny,
    uniform int num_markers
) {
    for (uniform int m = 0; m < num_markers; ++m) {
        uniform float mx = marker_x[m];
        uniform float my = marker_y[m];
        uniform float fx = marker_fx[m] * marker_ds[m];
        uniform float fy = marker_fy[m] * marker_ds[m];

        uniform int min_x = (int)floor(mx) - 2;
        uniform int max_x = (int)floor(mx) + 2;
        uniform int min_y = (int)floor(my) - 2;
        uniform int max_y = (int)floor(my) + 2;

        for (uniform int py = min_y; py <= max_y; ++py) {
            if (py >= 0 && py < ny) {
                uniform float wy = 0.25f * (1.0f + cos(0.5f * 3.14159265358979323846f * abs((uniform float)py - my)));
                uniform int row_offset = py * nx;

                foreach (px_idx = min_x ... max_x + 1) {
                    varying int px_wrapped = (px_idx % nx + nx) % nx;
                    varying float wx = peskin_kernel_1d((float)px_idx - mx);
                    varying float weight = wx * wy;

                    varying int idx = row_offset + px_wrapped;
                    euler_gx[idx] += weight * fx;
                    euler_gy[idx] += weight * fy;
                }
            }
        }
    }
}
