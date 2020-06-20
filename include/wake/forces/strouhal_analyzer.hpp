#pragma once

#include "wake/core/types.hpp"
#include <vector>
#include <complex>

namespace wake {

// Fast Fourier Transform (FFT) analyzer for periodic vortex shedding and Strouhal number extraction
class StrouhalAnalyzer {
public:
    StrouhalAnalyzer() = default;

    // Analyzes a lift coefficient time-series signal sampled with interval dt
    // Returns dominant vortex shedding frequency f_shed
    real_t analyze_frequency(const std::vector<real_t>& signal, real_t dt);

    // Calculates dimensionless Strouhal number St = f_shed * D / U_inf
    real_t compute_strouhal_number(real_t f_shed, real_t diameter, real_t u_inf);

    // Power spectral density distribution
    const std::vector<real_t>& frequencies() const noexcept { return frequencies_; }
    const std::vector<real_t>& power_spectrum() const noexcept { return power_spectrum_; }

private:
    void compute_radix2_fft(std::vector<std::complex<real_t>>& x);

    std::vector<real_t> frequencies_;
    std::vector<real_t> power_spectrum_;
};

} // namespace wake
