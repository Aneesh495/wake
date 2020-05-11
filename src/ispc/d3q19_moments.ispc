#include "wake/ispc/lbm_d3q19.isph"

// Exported ISPC kernel to compute macroscopic density and velocity from 3D D3Q19 distributions
export void d3q19_compute_macroscopic_moments(
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

        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            varying float fi = f[q * count + i];
            local_rho += fi;
            jx += fi * (float)d3q19_cx[q];
            jy += fi * (float)d3q19_cy[q];
            jz += fi * (float)d3q19_cz[q];
        }

        varying float inv_rho = 1.0f / local_rho;
        rho[i] = local_rho;
        ux[i] = jx * inv_rho;
        uy[i] = jy * inv_rho;
        uz[i] = jz * inv_rho;
    }
}

// Exported ISPC kernel to compute macroscopic moments with Guo 3D body forcing
export void d3q19_compute_macroscopic_moments_forced(
    const uniform float f[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float uz[],
    uniform float gx,
    uniform float gy,
    uniform float gz,
    uniform int count
) {
    foreach (i = 0 ... count) {
        varying float local_rho = 0.0f;
        varying float jx = 0.0f;
        varying float jy = 0.0f;
        varying float jz = 0.0f;

        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            varying float fi = f[q * count + i];
            local_rho += fi;
            jx += fi * (float)d3q19_cx[q];
            jy += fi * (float)d3q19_cy[q];
            jz += fi * (float)d3q19_cz[q];
        }

        varying float inv_rho = 1.0f / local_rho;
        rho[i] = local_rho;
        ux[i] = (jx + 0.5f * local_rho * gx) * inv_rho;
        uy[i] = (jy + 0.5f * local_rho * gy) * inv_rho;
        uz[i] = (jz + 0.5f * local_rho * gz) * inv_rho;
    }
}

// Compute 3D velocity magnitude across the lattice
export void d3q19_compute_velocity_magnitude(
    const uniform float ux[],
    const uniform float uy[],
    const uniform float uz[],
    uniform float speed[],
    uniform int count
) {
    foreach (i = 0 ... count) {
        varying float vx = ux[i];
        varying float vy = uy[i];
        varying float vz = uz[i];
        speed[i] = sqrt(vx * vx + vy * vy + vz * vz);
    }
}
