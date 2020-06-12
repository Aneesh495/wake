#include "wake/ispc/lbm_d3q27.isph"

// Exported ISPC kernel to compute macroscopic density and velocity from D3Q27 distributions
export void d3q27_compute_macroscopic_moments(
    const uniform float f[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float uz[],
    uniform int count
) {
    foreach (i = 0 ... count) {
        varying float local_rho = 0.0f;
        varying float jx = 0.0f;
        varying float jy = 0.0f;
        varying float jz = 0.0f;

        for (uniform int q = 0; q < D3Q27_Q; ++q) {
            varying float fi = f[q * count + i];
            local_rho += fi;
            jx += fi * (float)d3q27_cx[q];
            jy += fi * (float)d3q27_cy[q];
            jz += fi * (float)d3q27_cz[q];
        }

        varying float inv_rho = 1.0f / local_rho;
        rho[i] = local_rho;
        ux[i] = jx * inv_rho;
        uy[i] = jy * inv_rho;
        uz[i] = jz * inv_rho;
    }
}
