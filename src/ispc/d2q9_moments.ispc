#include "wake/ispc/lbm_d2q9.isph"

// Exported ISPC kernel to compute macroscopic density and velocity from D2Q9 distributions
// Distribution memory layout: Structure of Arrays (SoA) where f is arranged as Q blocks of size count.
export void d2q9_compute_macroscopic_moments(
    const uniform float f[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform int count
) {
    foreach (i = 0 ... count) {
        // Gather distribution values for all 9 discrete velocities
        varying float f0 = f[0 * count + i];
        varying float f1 = f[1 * count + i];
        varying float f2 = f[2 * count + i];
        varying float f3 = f[3 * count + i];
        varying float f4 = f[4 * count + i];
        varying float f5 = f[5 * count + i];
        varying float f6 = f[6 * count + i];
        varying float f7 = f[7 * count + i];
        varying float f8 = f[8 * count + i];

        // Zero-th moment: local macroscopic density
        varying float local_rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
        varying float inv_rho = 1.0f / local_rho;

        // First moment: momentum vector (rho * u)
        varying float jx = (f1 - f3) + (f5 - f6 - f7 + f8);
        varying float jy = (f2 - f4) + (f5 + f6 - f7 - f8);

        rho[i] = local_rho;
        ux[i] = jx * inv_rho;
        uy[i] = jy * inv_rho;
    }
}

// Compute velocity magnitude across the lattice
export void d2q9_compute_velocity_magnitude(
    const uniform float ux[],
    const uniform float uy[],
    uniform float speed[],
    uniform int count
) {
    foreach (i = 0 ... count) {
        varying float vx = ux[i];
        varying float vy = uy[i];
        speed[i] = sqrt(vx * vx + vy * vy);
    }
}
