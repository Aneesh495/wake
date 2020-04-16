#include "wake/ispc/lbm_d2q9.isph"

// Exported ISPC kernel for D2Q9 Multiple-Relaxation-Time (MRT) collision operator
// Transforms distribution functions into orthogonal Hermite moment space (Lallemand & Luo),
// relaxes moments independently according to physical dissipation rates, and transforms back.
export void d2q9_collide_mrt(
    const uniform float f_in[],
    uniform float f_out[],
    const uniform uint8 cell_type[],
    uniform float rho[],
    uniform float ux[],
    uniform float uy[],
    uniform float tau_shear,
    uniform float tau_bulk,
    uniform float gx,
    uniform float gy,
    uniform int count
) {
    uniform float s_rho = 0.0f;
    uniform float s_e   = 1.0f / tau_bulk;
    uniform float s_eps = 1.4f;
    uniform float s_jx  = 0.0f;
    uniform float s_q   = 1.6f;
    uniform float s_jy  = 0.0f;
    uniform float s_p   = 1.0f / tau_shear;

    foreach (i = 0 ... count) {
        if (cell_type[i] == CELL_SOLID) {
            for (uniform int q = 0; q < D2Q9_Q; ++q) {
                f_out[q * count + i] = f_in[q * count + i];
            }
            continue;
        }

        varying float f0 = f_in[0 * count + i];
        varying float f1 = f_in[1 * count + i];
        varying float f2 = f_in[2 * count + i];
        varying float f3 = f_in[3 * count + i];
        varying float f4 = f_in[4 * count + i];
        varying float f5 = f_in[5 * count + i];
        varying float f6 = f_in[6 * count + i];
        varying float f7 = f_in[7 * count + i];
        varying float f8 = f_in[8 * count + i];

        // 1. Forward orthogonal moment transformation: m = M * f
        varying float m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
        varying float m1 = -4.0f * f0 - (f1 + f2 + f3 + f4) + 2.0f * (f5 + f6 + f7 + f8);
        varying float m2 =  4.0f * f0 - 2.0f * (f1 + f2 + f3 + f4) + (f5 + f6 + f7 + f8);
        varying float m3 = (f1 - f3) + (f5 - f6 - f7 + f8);
        varying float m4 = -2.0f * (f1 - f3) + (f5 - f6 - f7 + f8);
        varying float m5 = (f2 - f4) + (f5 + f6 - f7 - f8);
        varying float m6 = -2.0f * (f2 - f4) + (f5 + f6 - f7 - f8);
        varying float m7 = (f1 - f2 + f3 - f4);
        varying float m8 = (f5 - f6 + f7 - f8);

        varying float local_rho = m0;
        varying float inv_rho = 1.0f / local_rho;

        // Guo velocity correction from body force
        varying float local_ux = (m3 + 0.5f * local_rho * gx) * inv_rho;
        varying float local_uy = (m5 + 0.5f * local_rho * gy) * inv_rho;

        rho[i] = local_rho;
        ux[i] = local_ux;
        uy[i] = local_uy;

        varying float ux2 = local_ux * local_ux;
        varying float uy2 = local_uy * local_uy;
        varying float u2  = ux2 + uy2;

        // 2. Equilibrium moments
        varying float meq0 = local_rho;
        varying float meq1 = local_rho * (-2.0f + 3.0f * u2);
        varying float meq2 = local_rho * (1.0f - 3.0f * u2);
        varying float meq3 = local_rho * local_ux;
        varying float meq4 = -local_rho * local_ux;
        varying float meq5 = local_rho * local_uy;
        varying float meq6 = -local_rho * local_uy;
        varying float meq7 = local_rho * (ux2 - uy2);
        varying float meq8 = local_rho * (local_ux * local_uy);

        // 3. Forcing terms in moment space
        varying float F_pre = (1.0f - 0.5f * s_p);
        varying float F_m3 = local_rho * gx;
        varying float F_m5 = local_rho * gy;
        varying float F_m1 = 6.0f * (local_ux * F_m3 + local_uy * F_m5);
        varying float F_m2 = -6.0f * (local_ux * F_m3 + local_uy * F_m5);
        varying float F_m4 = -F_m3;
        varying float F_m6 = -F_m5;
        varying float F_m7 = 2.0f * (local_ux * F_m3 - local_uy * F_m5);
        varying float F_m8 = (local_ux * F_m5 + local_uy * F_m3);

        // 4. Relax moments
        m0 = m0 - s_rho * (m0 - meq0);
        m1 = m1 - s_e   * (m1 - meq1) + F_pre * F_m1;
        m2 = m2 - s_eps * (m2 - meq2) + F_pre * F_m2;
        m3 = m3 - s_jx  * (m3 - meq3) + F_m3;
        m4 = m4 - s_q   * (m4 - meq4) + F_pre * F_m4;
        m5 = m5 - s_jy  * (m5 - meq5) + F_m5;
        m6 = m6 - s_q   * (m6 - meq6) + F_pre * F_m6;
        m7 = m7 - s_p   * (m7 - meq7) + F_pre * F_m7;
        m8 = m8 - s_p   * (m8 - meq8) + F_pre * F_m8;

        // 5. Inverse transformation: f = M^(-1) * m
        // Normalization weights:
        // m0: 1/9, m1: 1/36, m2: 1/36, m3: 1/6, m4: 1/12, m5: 1/6, m6: 1/12, m7: 1/4, m8: 1/4
        f_out[0 * count + i] = (1.0f / 9.0f) * m0 - (1.0f / 9.0f) * m1 + (1.0f / 9.0f) * m2;

        f_out[1 * count + i] = (1.0f / 9.0f) * m0 - (1.0f / 36.0f) * m1 - (1.0f / 18.0f) * m2
                             + (1.0f / 6.0f) * m3 - (1.0f / 6.0f) * m4 + (0.25f) * m7;

        f_out[2 * count + i] = (1.0f / 9.0f) * m0 - (1.0f / 36.0f) * m1 - (1.0f / 18.0f) * m2
                             + (1.0f / 6.0f) * m5 - (1.0f / 6.0f) * m6 - (0.25f) * m7;

        f_out[3 * count + i] = (1.0f / 9.0f) * m0 - (1.0f / 36.0f) * m1 - (1.0f / 18.0f) * m2
                             - (1.0f / 6.0f) * m3 + (1.0f / 6.0f) * m4 + (0.25f) * m7;

        f_out[4 * count + i] = (1.0f / 9.0f) * m0 - (1.0f / 36.0f) * m1 - (1.0f / 18.0f) * m2
                             - (1.0f / 6.0f) * m5 + (1.0f / 6.0f) * m6 - (0.25f) * m7;

        f_out[5 * count + i] = (1.0f / 9.0f) * m0 + (1.0f / 18.0f) * m1 + (1.0f / 36.0f) * m2
                             + (1.0f / 6.0f) * m3 + (1.0f / 12.0f) * m4
                             + (1.0f / 6.0f) * m5 + (1.0f / 12.0f) * m6 + (0.25f) * m8;

        f_out[6 * count + i] = (1.0f / 9.0f) * m0 + (1.0f / 18.0f) * m1 + (1.0f / 36.0f) * m2
                             - (1.0f / 6.0f) * m3 - (1.0f / 12.0f) * m4
                             + (1.0f / 6.0f) * m5 + (1.0f / 12.0f) * m6 - (0.25f) * m8;

        f_out[7 * count + i] = (1.0f / 9.0f) * m0 + (1.0f / 18.0f) * m1 + (1.0f / 36.0f) * m2
                             - (1.0f / 6.0f) * m3 - (1.0f / 12.0f) * m4
                             - (1.0f / 6.0f) * m5 - (1.0f / 12.0f) * m6 + (0.25f) * m8;

        f_out[8 * count + i] = (1.0f / 9.0f) * m0 + (1.0f / 18.0f) * m1 + (1.0f / 36.0f) * m2
                             + (1.0f / 6.0f) * m3 + (1.0f / 12.0f) * m4
                             - (1.0f / 6.0f) * m5 - (1.0f / 12.0f) * m6 - (0.25f) * m8;
    }
}
