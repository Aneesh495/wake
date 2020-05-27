#include "wake/ispc/lbm_d3q19.isph"

// Exported ISPC kernel for 3D D3Q19 Smagorinsky Large-Eddy Simulation (LES)
// Evaluates the full 3D non-equilibrium momentum flux tensor:
// Pi_xx, Pi_yy, Pi_zz, Pi_xy, Pi_yz, Pi_xz
// and calculates local effective turbulent relaxation time tau_eff
export void d3q19_compute_smagorinsky_tau(
    const uniform float f[],
    const uniform float rho[],
    const uniform float ux[],
    const uniform float uy[],
    const uniform float uz[],
    uniform float tau_eff[],
    uniform float tau_0,
    uniform float smagorinsky_constant,
    uniform int count
) {
    uniform float c_smag_sq = smagorinsky_constant * smagorinsky_constant;
    uniform float prefactor = 18.0f * 1.4142135623730951f * c_smag_sq;

    foreach (i = 0 ... count) {
        varying float local_rho = rho[i];
        varying float local_ux = ux[i];
        varying float local_uy = uy[i];
        varying float local_uz = uz[i];
        varying float u_sq = local_ux * local_ux + local_uy * local_uy + local_uz * local_uz;

        varying float pi_xx = 0.0f;
        varying float pi_yy = 0.0f;
        varying float pi_zz = 0.0f;
        varying float pi_xy = 0.0f;
        varying float pi_yz = 0.0f;
        varying float pi_xz = 0.0f;

        for (uniform int q = 0; q < D3Q19_Q; ++q) {
            varying float fi = f[q * count + i];

            varying float ci_x = (float)d3q19_cx[q];
            varying float ci_y = (float)d3q19_cy[q];
            varying float ci_z = (float)d3q19_cz[q];
            varying float ci_u = ci_x * local_ux + ci_y * local_uy + ci_z * local_uz;

            varying float feq = d3q19_w[q] * local_rho * (
                1.0f + LBM_INV_CS2 * ci_u + LBM_INV_2_CS4 * ci_u * ci_u - LBM_INV_2_CS2 * u_sq
            );

            varying float f_neq = fi - feq;

            pi_xx += ci_x * ci_x * f_neq;
            pi_yy += ci_y * ci_y * f_neq;
            pi_zz += ci_z * ci_z * f_neq;
            pi_xy += ci_x * ci_y * f_neq;
            pi_yz += ci_y * ci_z * f_neq;
            pi_xz += ci_x * ci_z * f_neq;
        }

        // 3D Frobenius norm: Q = sqrt(2 * (Pi_xx^2 + Pi_yy^2 + Pi_zz^2 + 2 * (Pi_xy^2 + Pi_yz^2 + Pi_xz^2)))
        varying float tensor_sum = (pi_xx * pi_xx + pi_yy * pi_yy + pi_zz * pi_zz)
                                 + 2.0f * (pi_xy * pi_xy + pi_yz * pi_yz + pi_xz * pi_xz);
        varying float Q_stress = sqrt(2.0f * tensor_sum);

        // Solve quadratic equation: tau_eff = (tau_0 + sqrt(tau_0^2 + 18 * sqrt(2) * C_s^2 / rho * Q)) / 2
        varying float disc = tau_0 * tau_0 + (prefactor / local_rho) * Q_stress;
        tau_eff[i] = 0.5f * (tau_0 + sqrt(disc));
    }
}
