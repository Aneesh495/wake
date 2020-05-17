#include "wake/ispc/lbm_d3q19.isph"

// Exported ISPC kernel for 3D D3Q19 Multiple-Relaxation-Time (MRT) collision operator
// Transforms populations into 19 orthogonal Hermite moments (d'Humieres et al.),
// relaxes each hydrodynamic and kinetic ghost mode with its physical relaxation rate,
// and projects back to distribution space.
export void d3q19_collide_mrt(
    const uniform float f_in[],
    uniform float f_out[],
    const uniform uint8 cell_type[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float uz[],
    uniform float tau_shear,
    uniform float tau_bulk,
    uniform float gx,
    uniform float gy,
    uniform float gz,
    uniform int count
) {
    uniform float s_rho = 0.0f;
    uniform float s_e   = 1.0f / tau_bulk; // Energy mode governs bulk viscosity
    uniform float s_eps = 1.4f;
    uniform float s_j   = 0.0f;            // Momentum conservation
    uniform float s_q   = 1.2f;            // Heat flux modes
    uniform float s_nu  = 1.0f / tau_shear;// Viscous shear stress tensor modes
    uniform float s_pi  = 1.4f;
    uniform float s_m   = 1.98f;

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID) {
            for (uniform int q = 0; q < D3Q19_Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        varying float f[D3Q19_Q];
        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            f[q] = f_in[q * count + i];
        }

        // 1. Forward 19-moment transformation: m = M * f
        // m0: Density
        varying float m0 = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6]
                         + f[7] + f[8] + f[9] + f[10] + f[11] + f[12]
                         + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];

        // m1: Energy (e)
        varying float m1 = -30.0f * f[0]
            - 11.0f * (f[1] + f[2] + f[3] + f[4] + f[5] + f[6])
            + 8.0f * (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]);

        // m2: Energy squared (epsilon)
        varying float m2 = 12.0f * f[0]
            - 4.0f * (f[1] + f[2] + f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]);

        // m3, m5, m7: Momentum components (jx, jy, jz)
        varying float m3 = (f[1] - f[2]) + (f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14]);
        varying float m5 = (f[3] - f[4]) + (f[7] + f[8] - f[9] - f[10] + f[15] - f[16] + f[17] - f[18]);
        varying float m7 = (f[5] - f[6]) + (f[11] + f[12] - f[13] - f[14] + f[15] + f[16] - f[17] - f[18]);

        // m4, m6, m8: Heat flux components (qx, qy, qz)
        varying float m4 = -4.0f * (f[1] - f[2]) + (f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14]);
        varying float m6 = -4.0f * (f[3] - f[4]) + (f[7] + f[8] - f[9] - f[10] + f[15] - f[16] + f[17] - f[18]);
        varying float m8 = -4.0f * (f[5] - f[6]) + (f[11] + f[12] - f[13] - f[14] + f[15] + f[16] - f[17] - f[18]);

        // m9, m11: Diagonal stress components
        varying float m9 = 2.0f * (f[1] + f[2]) - (f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14]) - 2.0f * (f[15] + f[16] + f[17] + f[18]);

        varying float m11 = (f[3] + f[4]) - (f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10]) - (f[11] + f[12] + f[13] + f[14]);

        // m10, m12: High-order stress components
        varying float m10 = -4.0f * (f[1] + f[2]) + 2.0f * (f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14]) - 2.0f * (f[15] + f[16] + f[17] + f[18]);

        varying float m12 = -2.0f * (f[3] + f[4]) + 2.0f * (f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10]) - (f[11] + f[12] + f[13] + f[14]);

        // m13, m14, m15: Off-diagonal shear stress components (pxy, pyz, pxz)
        varying float m13 = (f[7] - f[8] - f[9] + f[10]);
        varying float m14 = (f[15] - f[16] - f[17] + f[18]);
        varying float m15 = (f[11] - f[12] - f[13] + f[14]);

        // m16, m17, m18: Third-order antisymmetric modes
        varying float m16 = (f[7] - f[8] + f[9] - f[10] - f[11] + f[12] - f[13] + f[14]);
        varying float m17 = (-f[7] - f[8] + f[9] + f[10] + f[15] - f[16] + f[17] - f[18]);
        varying float m18 = (-f[11] - f[12] + f[13] + f[14] + f[15] + f[16] - f[17] - f[18]);

        varying float local_rho = m0;
        varying float inv_rho = 1.0f / local_rho;

        // Guo velocity correction from body force
        varying float local_ux = (m3 + 0.5f * local_rho * gx) * inv_rho;
        varying float local_uy = (m5 + 0.5f * local_rho * gy) * inv_rho;
        varying float local_uz = (m7 + 0.5f * local_rho * gz) * inv_rho;

        rho[i] = local_rho;
        ux[i] = local_ux;
        uy[i] = local_uy;
        uz[i] = local_uz;

        varying float ux2 = local_ux * local_ux;
        varying float uy2 = local_uy * local_uy;
        varying float uz2 = local_uz * local_uz;
        varying float u2  = ux2 + uy2 + uz2;

        // 2. Equilibrium moments
        varying float meq0 = local_rho;
        varying float meq1 = local_rho * (-11.0f + 19.0f * u2);
        varying float meq2 = local_rho * (3.0f - 11.0f / 2.0f * u2);
        varying float meq3 = local_rho * local_ux;
        varying float meq4 = -2.0f / 3.0f * local_rho * local_ux;
        varying float meq5 = local_rho * local_uy;
        varying float meq6 = -2.0f / 3.0f * local_rho * local_uy;
        varying float meq7 = local_rho * local_uz;
        varying float meq8 = -2.0f / 3.0f * local_rho * local_uz;
        varying float meq9 = local_rho * (2.0f * ux2 - uy2 - uz2);
        varying float meq10 = -0.5f * meq9;
        varying float meq11 = local_rho * (uy2 - uz2);
        varying float meq12 = -0.5f * meq11;
        varying float meq13 = local_rho * (local_ux * local_uy);
        varying float meq14 = local_rho * (local_uy * local_uz);
        varying float meq15 = local_rho * (local_ux * local_uz);
        varying float meq16 = 0.0f;
        varying float meq17 = 0.0f;
        varying float meq18 = 0.0f;

        // 3. Relax moments
        m0 = m0 - s_rho * (m0 - meq0);
        m1 = m1 - s_e   * (m1 - meq1);
        m2 = m2 - s_eps * (m2 - meq2);
        m3 = m3 - s_j   * (m3 - meq3) + local_rho * gx;
        m4 = m4 - s_q   * (m4 - meq4) - (2.0f / 3.0f) * local_rho * gx;
        m5 = m5 - s_j   * (m5 - meq5) + local_rho * gy;
        m6 = m6 - s_q   * (m6 - meq6) - (2.0f / 3.0f) * local_rho * gy;
        m7 = m7 - s_j   * (m7 - meq7) + local_rho * gz;
        m8 = m8 - s_q   * (m8 - meq8) - (2.0f / 3.0f) * local_rho * gz;
        m9  = m9  - s_nu * (m9  - meq9);
        m10 = m10 - s_pi * (m10 - meq10);
        m11 = m11 - s_nu * (m11 - meq11);
        m12 = m12 - s_pi * (m12 - meq12);
        m13 = m13 - s_nu * (m13 - meq13);
        m14 = m14 - s_nu * (m14 - meq14);
        m15 = m15 - s_nu * (m15 - meq15);
        m16 = m16 - s_m  * (m16 - meq16);
        m17 = m17 - s_m  * (m17 - meq17);
        m18 = m18 - s_m  * (m18 - meq18);

        // 4. Inverse transformation: f = M^(-1) * m
        // Normalization denominators of M rows:
        // 19, 2394, 252, 10, 40, 10, 40, 10, 40, 162, 252, 18, 72, 4, 4, 4, 40, 40, 40
        f_out[0 * count + i] = (1.0f / 19.0f) * m0 - (5.0f / 399.0f) * m1 + (1.0f / 21.0f) * m2;

        f_out[1 * count + i] = (1.0f / 19.0f) * m0 - (11.0f / 2394.0f) * m1 - (1.0f / 63.0f) * m2
                             + (0.1f) * m3 - (0.1f) * m4 + (1.0f / 81.0f) * m9 - (1.0f / 63.0f) * m10;

        f_out[2 * count + i] = (1.0f / 19.0f) * m0 - (11.0f / 2394.0f) * m1 - (1.0f / 63.0f) * m2
                             - (0.1f) * m3 + (0.1f) * m4 + (1.0f / 81.0f) * m9 - (1.0f / 63.0f) * m10;

        f_out[3 * count + i] = (1.0f / 19.0f) * m0 - (11.0f / 2394.0f) * m1 - (1.0f / 63.0f) * m2
                             + (0.1f) * m5 - (0.1f) * m6 - (1.0f / 162.0f) * m9 + (1.0f / 126.0f) * m10
                             + (1.0f / 18.0f) * m11 - (1.0f / 36.0f) * m12;

        f_out[4 * count + i] = (1.0f / 19.0f) * m0 - (11.0f / 2394.0f) * m1 - (1.0f / 63.0f) * m2
                             - (0.1f) * m5 + (0.1f) * m6 - (1.0f / 162.0f) * m9 + (1.0f / 126.0f) * m10
                             + (1.0f / 18.0f) * m11 - (1.0f / 36.0f) * m12;

        f_out[5 * count + i] = (1.0f / 19.0f) * m0 - (11.0f / 2394.0f) * m1 - (1.0f / 63.0f) * m2
                             + (0.1f) * m7 - (0.1f) * m8 - (1.0f / 162.0f) * m9 + (1.0f / 126.0f) * m10
                             - (1.0f / 18.0f) * m11 + (1.0f / 36.0f) * m12;

        f_out[6 * count + i] = (1.0f / 19.0f) * m0 - (11.0f / 2394.0f) * m1 - (1.0f / 63.0f) * m2
                             - (0.1f) * m7 + (0.1f) * m8 - (1.0f / 162.0f) * m9 + (1.0f / 126.0f) * m10
                             - (1.0f / 18.0f) * m11 + (1.0f / 36.0f) * m12;

        f_out[7 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                             + (1.0f / 20.0f) * (m3 + m5) + (1.0f / 40.0f) * (m4 + m6)
                             + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                             + (1.0f / 36.0f) * m11 + (1.0f / 72.0f) * m12 + (0.25f) * m13;

        f_out[8 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                             + (1.0f / 20.0f) * (-m3 + m5) + (1.0f / 40.0f) * (-m4 + m6)
                             + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                             + (1.0f / 36.0f) * m11 + (1.0f / 72.0f) * m12 - (0.25f) * m13;

        f_out[9 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                             + (1.0f / 20.0f) * (m3 - m5) + (1.0f / 40.0f) * (m4 - m6)
                             + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                             + (1.0f / 36.0f) * m11 + (1.0f / 72.0f) * m12 - (0.25f) * m13;

        f_out[10 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (-m3 - m5) + (1.0f / 40.0f) * (-m4 - m6)
                              + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                              + (1.0f / 36.0f) * m11 + (1.0f / 72.0f) * m12 + (0.25f) * m13;

        f_out[11 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (m3 + m7) + (1.0f / 40.0f) * (m4 + m8)
                              + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                              - (1.0f / 36.0f) * m11 - (1.0f / 72.0f) * m12 + (0.25f) * m15;

        f_out[12 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (-m3 + m7) + (1.0f / 40.0f) * (-m4 + m8)
                              + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                              - (1.0f / 36.0f) * m11 - (1.0f / 72.0f) * m12 - (0.25f) * m15;

        f_out[13 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (m3 - m7) + (1.0f / 40.0f) * (m4 - m8)
                              + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                              - (1.0f / 36.0f) * m11 - (1.0f / 72.0f) * m12 - (0.25f) * m15;

        f_out[14 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (-m3 - m7) + (1.0f / 40.0f) * (-m4 - m8)
                              + (1.0f / 162.0f) * m9 + (1.0f / 252.0f) * m10
                              - (1.0f / 36.0f) * m11 - (1.0f / 72.0f) * m12 + (0.25f) * m15;

        f_out[15 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (m5 + m7) + (1.0f / 40.0f) * (m6 + m8)
                              - (1.0f / 81.0f) * m9 - (1.0f / 126.0f) * m10 + (0.25f) * m14;

        f_out[16 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (-m5 + m7) + (1.0f / 40.0f) * (-m6 + m8)
                              - (1.0f / 81.0f) * m9 - (1.0f / 126.0f) * m10 - (0.25f) * m14;

        f_out[17 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (m5 - m7) + (1.0f / 40.0f) * (m6 - m8)
                              - (1.0f / 81.0f) * m9 - (1.0f / 126.0f) * m10 - (0.25f) * m14;

        f_out[18 * count + i] = (1.0f / 19.0f) * m0 + (4.0f / 1197.0f) * m1 + (1.0f / 252.0f) * m2
                              + (1.0f / 20.0f) * (-m5 - m7) + (1.0f / 40.0f) * (-m6 - m8)
                              - (1.0f / 81.0f) * m9 - (1.0f / 126.0f) * m10 + (0.25f) * m14;
    }
}
