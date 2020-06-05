#include "wake/ispc/lbm_d3q19.isph"

// Exact analytical inverse transformation matrix Minv = inv(M) for d'Humieres D3Q19 MRT
static const uniform float d3q19_Minv[D3Q19_Q][D3Q19_Q] = {
    {0.0526316f, -0.0125313f, 0.0476190f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, -0.0045948f, -0.0158730f, 0.1000000f, -0.1000000f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0555556f, -0.0555556f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, -0.0045948f, -0.0158730f, -0.1000000f, 0.1000000f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0555556f, -0.0555556f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, -0.0045948f, -0.0158730f, 0.0f, 0.0f, 0.1000000f, -0.1000000f, 0.0f, 0.0f, -0.0277778f, 0.0277778f, 0.0833333f, -0.0833333f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, -0.0045948f, -0.0158730f, 0.0f, 0.0f, -0.1000000f, 0.1000000f, 0.0f, 0.0f, -0.0277778f, 0.0277778f, 0.0833333f, -0.0833333f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, -0.0045948f, -0.0158730f, 0.0f, 0.0f, 0.0f, 0.0f, 0.1000000f, -0.1000000f, -0.0277778f, 0.0277778f, -0.0833333f, 0.0833333f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, -0.0045948f, -0.0158730f, 0.0f, 0.0f, 0.0f, 0.0f, -0.1000000f, 0.1000000f, -0.0277778f, 0.0277778f, -0.0833333f, 0.0833333f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.1000000f, 0.0250000f, 0.1000000f, 0.0250000f, 0.0f, 0.0f, 0.0277778f, 0.0138889f, 0.0833333f, 0.0416667f, 0.2500000f, 0.0f, 0.0f, 0.1250000f, -0.1250000f, 0.0f},
    {0.0526316f, 0.0033417f, 0.0039683f, -0.1000000f, -0.0250000f, 0.1000000f, 0.0250000f, 0.0f, 0.0f, 0.0277778f, 0.0138889f, 0.0833333f, 0.0416667f, -0.2500000f, 0.0f, 0.0f, -0.1250000f, -0.1250000f, 0.0f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.1000000f, 0.0250000f, -0.1000000f, -0.0250000f, 0.0f, 0.0f, 0.0277778f, 0.0138889f, 0.0833333f, 0.0416667f, -0.2500000f, 0.0f, 0.0f, 0.1250000f, 0.1250000f, 0.0f},
    {0.0526316f, 0.0033417f, 0.0039683f, -0.1000000f, -0.0250000f, -0.1000000f, -0.0250000f, 0.0f, 0.0f, 0.0277778f, 0.0138889f, 0.0833333f, 0.0416667f, 0.2500000f, 0.0f, 0.0f, -0.1250000f, 0.1250000f, 0.0f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.1000000f, 0.0250000f, 0.0f, 0.0f, 0.1000000f, 0.0250000f, 0.0277778f, 0.0138889f, -0.0833333f, -0.0416667f, 0.0f, 0.0f, 0.2500000f, -0.1250000f, 0.0f, -0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, -0.1000000f, -0.0250000f, 0.0f, 0.0f, 0.1000000f, 0.0250000f, 0.0277778f, 0.0138889f, -0.0833333f, -0.0416667f, 0.0f, 0.0f, -0.2500000f, 0.1250000f, 0.0f, -0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.1000000f, 0.0250000f, 0.0f, 0.0f, -0.1000000f, -0.0250000f, 0.0277778f, 0.0138889f, -0.0833333f, -0.0416667f, 0.0f, 0.0f, -0.2500000f, -0.1250000f, 0.0f, 0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, -0.1000000f, -0.0250000f, 0.0f, 0.0f, -0.1000000f, -0.0250000f, 0.0277778f, 0.0138889f, -0.0833333f, -0.0416667f, 0.0f, 0.0f, 0.2500000f, 0.1250000f, 0.0f, 0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.0f, 0.0f, 0.1000000f, 0.0250000f, 0.1000000f, 0.0250000f, -0.0555556f, -0.0277778f, 0.0f, 0.0f, 0.0f, 0.2500000f, 0.0f, 0.0f, 0.1250000f, 0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.0f, 0.0f, -0.1000000f, -0.0250000f, 0.1000000f, 0.0250000f, -0.0555556f, -0.0277778f, 0.0f, 0.0f, 0.0f, -0.2500000f, 0.0f, 0.0f, -0.1250000f, 0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.0f, 0.0f, 0.1000000f, 0.0250000f, -0.1000000f, -0.0250000f, -0.0555556f, -0.0277778f, 0.0f, 0.0f, 0.0f, -0.2500000f, 0.0f, 0.0f, 0.1250000f, -0.1250000f},
    {0.0526316f, 0.0033417f, 0.0039683f, 0.0f, 0.0f, -0.1000000f, -0.0250000f, -0.1000000f, -0.0250000f, -0.0555556f, -0.0277778f, 0.0f, 0.0f, 0.0f, 0.2500000f, 0.0f, 0.0f, -0.1250000f, -0.1250000f}
};

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
    uniform float s_e   = 1.0f / tau_bulk;
    uniform float s_eps = 1.4f;
    uniform float s_j   = 0.0f;
    uniform float s_q   = 1.2f;
    uniform float s_nu  = 1.0f / tau_shear;
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

        // 1. Forward orthogonal moment transformation
        varying float m[D3Q19_Q];

        m[0] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6]
             + f[7] + f[8] + f[9] + f[10] + f[11] + f[12]
             + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];

        m[1] = -30.0f * f[0]
            - 11.0f * (f[1] + f[2] + f[3] + f[4] + f[5] + f[6])
            + 8.0f * (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]);

        m[2] = 12.0f * f[0]
            - 4.0f * (f[1] + f[2] + f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18]);

        m[3] = (f[1] - f[2]) + (f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14]);
        m[4] = -4.0f * (f[1] - f[2]) + (f[7] - f[8] + f[9] - f[10] + f[11] - f[12] + f[13] - f[14]);

        m[5] = (f[3] - f[4]) + (f[7] + f[8] - f[9] - f[10] + f[15] - f[16] + f[17] - f[18]);
        m[6] = -4.0f * (f[3] - f[4]) + (f[7] + f[8] - f[9] - f[10] + f[15] - f[16] + f[17] - f[18]);

        m[7] = (f[5] - f[6]) + (f[11] + f[12] - f[13] - f[14] + f[15] + f[16] - f[17] - f[18]);
        m[8] = -4.0f * (f[5] - f[6]) + (f[11] + f[12] - f[13] - f[14] + f[15] + f[16] - f[17] - f[18]);

        m[9] = 2.0f * (f[1] + f[2]) - (f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14]) - 2.0f * (f[15] + f[16] + f[17] + f[18]);

        m[10] = -4.0f * (f[1] + f[2]) + 2.0f * (f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14]) - 2.0f * (f[15] + f[16] + f[17] + f[18]);

        m[11] = (f[3] + f[4]) - (f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10]) - (f[11] + f[12] + f[13] + f[14]);

        m[12] = -2.0f * (f[3] + f[4]) + 2.0f * (f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10]) - (f[11] + f[12] + f[13] + f[14]);

        m[13] = (f[7] - f[8] - f[9] + f[10]);
        m[14] = (f[15] - f[16] - f[17] + f[18]);
        m[15] = (f[11] - f[12] - f[13] + f[14]);

        m[16] = (f[7] - f[8] + f[9] - f[10] - f[11] + f[12] - f[13] + f[14]);
        m[17] = (-f[7] - f[8] + f[9] + f[10] + f[15] - f[16] + f[17] - f[18]);
        m[18] = (-f[11] - f[12] + f[13] + f[14] + f[15] + f[16] - f[17] - f[18]);

        varying float local_rho = m[0];
        varying float inv_rho = 1.0f / local_rho;

        varying float local_ux = (m[3] + 0.5f * local_rho * gx) * inv_rho;
        varying float local_uy = (m[5] + 0.5f * local_rho * gy) * inv_rho;
        varying float local_uz = (m[7] + 0.5f * local_rho * gz) * inv_rho;

        rho[i] = local_rho;
        ux[i] = local_ux;
        uy[i] = local_uy;
        uz[i] = local_uz;

        varying float ux2 = local_ux * local_ux;
        varying float uy2 = local_uy * local_uy;
        varying float uz2 = local_uz * local_uz;
        varying float u2  = ux2 + uy2 + uz2;

        // 2. Equilibrium moments
        varying float meq[D3Q19_Q];
        meq[0] = local_rho;
        meq[1] = local_rho * (-11.0f + 19.0f * u2);
        meq[2] = local_rho * (3.0f - 11.0f / 2.0f * u2);
        meq[3] = local_rho * local_ux;
        meq[4] = -2.0f / 3.0f * local_rho * local_ux;
        meq[5] = local_rho * local_uy;
        meq[6] = -2.0f / 3.0f * local_rho * local_uy;
        meq[7] = local_rho * local_uz;
        meq[8] = -2.0f / 3.0f * local_rho * local_uz;
        meq[9] = local_rho * (2.0f * ux2 - uy2 - uz2);
        meq[10] = -0.5f * meq[9];
        meq[11] = local_rho * (uy2 - uz2);
        meq[12] = -0.5f * meq[11];
        meq[13] = local_rho * (local_ux * local_uy);
        meq[14] = local_rho * (local_uy * local_uz);
        meq[15] = local_rho * (local_ux * local_uz);
        meq[16] = 0.0f;
        meq[17] = 0.0f;
        meq[18] = 0.0f;

        // 3. Relax moments
        m[0] = m[0] - s_rho * (m[0] - meq[0]);
        m[1] = m[1] - s_e   * (m[1] - meq[1]);
        m[2] = m[2] - s_eps * (m[2] - meq[2]);
        m[3] = m[3] - s_j   * (m[3] - meq[3]) + local_rho * gx;
        m[4] = m[4] - s_q   * (m[4] - meq[4]) - (2.0f / 3.0f) * local_rho * gx;
        m[5] = m[5] - s_j   * (m[5] - meq[5]) + local_rho * gy;
        m[6] = m[6] - s_q   * (m[6] - meq[6]) - (2.0f / 3.0f) * local_rho * gy;
        m[7] = m[7] - s_j   * (m[7] - meq[7]) + local_rho * gz;
        m[8] = m[8] - s_q   * (m[8] - meq[8]) - (2.0f / 3.0f) * local_rho * gz;
        m[9]  = m[9]  - s_nu * (m[9]  - meq[9]);
        m[10] = m[10] - s_pi * (m[10] - meq[10]);
        m[11] = m[11] - s_nu * (m[11] - meq[11]);
        m[12] = m[12] - s_pi * (m[12] - meq[12]);
        m[13] = m[13] - s_nu * (m[13] - meq[13]);
        m[14] = m[14] - s_nu * (m[14] - meq[14]);
        m[15] = m[15] - s_nu * (m[15] - meq[15]);
        m[16] = m[16] - s_m  * (m[16] - meq[16]);
        m[17] = m[17] - s_m  * (m[17] - meq[17]);
        m[18] = m[18] - s_m  * (m[18] - meq[18]);

        // 4. Exact inverse transformation: f = Minv * m
        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            varying float val = 0.0f;
            for (uniform int k = 0; k < D3Q19_Q; ++k) {
                val += d3q19_Minv[q][k] * m[k];
            }
            f_out[q * count + i] = val;
        }
    }
}
