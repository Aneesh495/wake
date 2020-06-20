#include "wake/forces/strouhal_analyzer.hpp"
#include "wake/core/constants.hpp"
#include <cmath>
#include <algorithm>

namespace wake {

void StrouhalAnalyzer::compute_radix2_fft(std::vector<std::complex<real_t>>& a) {
    const std::size_t n = a.size();
    if (n <= 1) return;

    // Bit-reversal permutation
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // Cooley-Tukey butterfly stages
    for (std::size_t len = 2; len <= n; len <<= 1) {
        real_t angle = -TWO_PI / static_cast<real_t>(len);
        std::complex<real_t> wlen(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<real_t> w(1.0f, 0.0f);
            for (std::size_t j = 0; j < len / 2; ++j) {
                std::complex<real_t> u = a[i + j];
                std::complex<real_t> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

real_t StrouhalAnalyzer::analyze_frequency(const std::vector<real_t>& signal, real_t dt) {
    if (signal.size() < 16 || dt <= 0.0f) return 0.0f;

    // Remove mean DC offset
    real_t sum = 0.0f;
    for (real_t v : signal) sum += v;
    real_t mean = sum / static_cast<real_t>(signal.size());

    // Zero-pad to nearest power of 2
    std::size_t n = 1;
    while (n < signal.size()) n <<= 1;
    if (n > 4096) n = 4096; // Bound maximum FFT size

    std::vector<std::complex<real_t>> a(n, std::complex<real_t>(0.0f, 0.0f));
    std::size_t input_len = std::min(n, signal.size());

    // Apply Hann window to eliminate spectral leakage
    for (std::size_t i = 0; i < input_len; ++i) {
        real_t window = 0.5f * (1.0f - std::cos(TWO_PI * static_cast<real_t>(i) / static_cast<real_t>(input_len - 1)));
        a[i] = std::complex<real_t>((signal[i] - mean) * window, 0.0f);
    }

    compute_radix2_fft(a);

    // Compute power spectral density
    std::size_t half_n = n / 2;
    frequencies_.resize(half_n);
    power_spectrum_.resize(half_n);

    real_t df = 1.0f / (static_cast<real_t>(n) * dt);
    std::size_t peak_idx = 1;
    real_t max_power = 0.0f;

    // Ignore DC bin (idx = 0)
    for (std::size_t i = 1; i < half_n; ++i) {
        frequencies_[i] = static_cast<real_t>(i) * df;
        real_t p = std::norm(a[i]);
        power_spectrum_[i] = p;
        if (p > max_power) {
            max_power = p;
            peak_idx = i;
        }
    }

    return frequencies_[peak_idx];
}

real_t StrouhalAnalyzer::compute_strouhal_number(real_t f_shed, real_t diameter, real_t u_inf) {
    if (u_inf <= 1e-8f) return 0.0f;
    return (f_shed * diameter) / u_inf;
}

} // namespace wake
