#include "wake/ispc/lbm_d2q9.isph"

// Exported ISPC kernel for Smagorinsky Large-Eddy Simulation (LES) in D2Q9
// Computes local non-equilibrium stress tensor Pi_neq directly from distributions
// and solves the quadratic Smagorinsky relation for local effective relaxation time tau_eff
export void d2q9_compute_smagorinsky_tau(
    const uniform float f[],
    const uniform float rho[],
    const uniform float ux[],
    const uniform float uy[],
    uniform float tau_eff[],
    uniform float tau_0,
    uniform float smagorinsky_constant,
    uniform int count
) {
    // Constant C_smag^2 * Delta^2
    uniform float c_smag_sq = smagorinsky_constant * smagorinsky_constant;
    uniform float prefactor = 18.0f * 1.4142135623730951f * c_smag_sq;

    foreach (i = 0 ... count) {
        varying float local_rho = rho[i];
        varying float local_ux = ux[i];
        varying float local_uy = uy[i];
        varying float u_sq = local_ux * local_ux + local_uy * local_uy;

        // Components of non-equilibrium momentum flux tensor: Pi_xx, Pi_xy, Pi_yy
        varying float pi_xx = 0.0f;
        varying float pi_xy = 0.0f;
        varying float pi_yy = 0.0f;

        for (uniform int q = 0; q < D2Q9_Q; ++q) {
            varying float fi = f[q * count + i];

            varying float ci_x = (float)d2q9_cx[q];
            varying float ci_y = (float)d2q9_cy[q];
            varying float ci_u = ci_x * local_ux + ci_y * local_uy;

            varying float feq = d2q9_w[q] * local_rho * (
                1.0f + LBM_INV_CS2 * ci_u + LBM_INV_2_CS4 * ci_u * ci_u - LBM_INV_2_CS2 * u_sq
            );

            varying float f_neq = fi - feq;

            pi_xx += ci_x * ci_x * f_neq;
            pi_xy += ci_x * ci_y * f_neq;
            pi_yy += ci_y * ci_y * f_neq;
        }

        // Tensor Frobenius norm: Q = sqrt(2 * (Pi_xx^2 + 2 * Pi_xy^2 + Pi_yy^2))
        varying float tensor_sum = pi_xx * pi_xx + 2.0f * pi_xy * pi_xy + pi_yy * pi_yy;
        varying float Q_stress = sqrt(2.0f * tensor_sum);

        // Solve quadratic equation: tau_eff = (tau_0 + sqrt(tau_0^2 + 18 * sqrt(2) * C_s^2 / rho * Q)) / 2
        varying float disc = tau_0 * tau_0 + (prefactor / local_rho) * Q_stress;
        tau_eff[i] = 0.5f * (tau_0 + sqrt(disc));
    }
}
