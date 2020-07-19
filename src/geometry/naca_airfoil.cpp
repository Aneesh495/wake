#include "wake/geometry/naca_airfoil.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace wake {

NacaAirfoil::NacaAirfoil(
    const std::string& code,
    real_t chord,
    real_t angle_of_attack_deg,
    size_t num_points
) {
    configure(code, chord, angle_of_attack_deg, num_points);
}

void NacaAirfoil::configure(
    const std::string& code,
    real_t chord,
    real_t angle_of_attack_deg,
    size_t num_points
) {
    code_ = code;
    chord_ = chord;
    aoa_deg_ = angle_of_attack_deg;

    if (code.length() != 4) {
        throw std::invalid_argument("NACA code must be 4 digits, e.g. 0012");
    }

    int m_digit = code[0] - '0';
    int p_digit = code[1] - '0';
    int t_digits = std::stoi(code.substr(2, 2));

    camber_ = static_cast<real_t>(m_digit) / 100.0f;
    camber_pos_ = static_cast<real_t>(p_digit) / 10.0f;
    thickness_ = static_cast<real_t>(t_digits) / 100.0f;

    generate_coordinates(num_points);
}

void NacaAirfoil::generate_coordinates(size_t num_points) {
    polygon_.clear();
    polygon_.reserve(num_points * 2);

    real_t alpha_rad = -aoa_deg_ * (PI / 180.0f);
    real_t cos_a = std::cos(alpha_rad);
    real_t sin_a = std::sin(alpha_rad);

    std::vector<AirfoilPoint> upper;
    std::vector<AirfoilPoint> lower;
    upper.reserve(num_points);
    lower.reserve(num_points);

    // Cosine distribution of chordwise stations x in [0, 1]
    for (size_t i = 0; i <= num_points; ++i) {
        real_t beta = static_cast<real_t>(i) * PI / static_cast<real_t>(num_points);
        real_t x = 0.5f * (1.0f - std::cos(beta));

        // Thickness distribution
        real_t yt = 5.0f * thickness_ * (
            0.2969f * std::sqrt(x) -
            0.1260f * x -
            0.3516f * x * x +
            0.2843f * x * x * x -
            0.1036f * x * x * x * x // closed trailing edge
        );

        // Camber line and gradient
        real_t yc = 0.0f;
        real_t dyc_dx = 0.0f;

        if (camber_pos_ > 1e-6f) {
            if (x < camber_pos_) {
                yc = (camber_ / (camber_pos_ * camber_pos_)) * (2.0f * camber_pos_ * x - x * x);
                dyc_dx = (2.0f * camber_ / (camber_pos_ * camber_pos_)) * (camber_pos_ - x);
            } else {
                real_t one_minus_p = 1.0f - camber_pos_;
                yc = (camber_ / (one_minus_p * one_minus_p)) * ((1.0f - 2.0f * camber_pos_) + 2.0f * camber_pos_ * x - x * x);
                dyc_dx = (2.0f * camber_ / (one_minus_p * one_minus_p)) * (camber_pos_ - x);
            }
        }

        real_t theta = std::atan(dyc_dx);
        real_t cos_th = std::cos(theta);
        real_t sin_th = std::sin(theta);

        // Upper and lower profile points
        real_t xu = (x - yt * sin_th) * chord_;
        real_t yu = (yc + yt * cos_th) * chord_;

        real_t xl = (x + yt * sin_th) * chord_;
        real_t yl = (yc - yt * cos_th) * chord_;

        // Apply angle of attack rotation
        upper.push_back({xu * cos_a - yu * sin_a, xu * sin_a + yu * cos_a});
        lower.push_back({xl * cos_a - yl * sin_a, xl * sin_a + yl * cos_a});
    }

    // Assemble closed polygon: Leading Edge -> Upper surface -> Trailing Edge -> Lower surface -> Leading Edge
    // Upper from LE (i=0) to TE (i=num_points)
    for (size_t i = 0; i <= num_points; ++i) {
        polygon_.push_back(upper[i]);
    }
    // Lower from TE back to LE
    for (size_t i = num_points; i > 0; --i) {
        polygon_.push_back(lower[i]);
    }
}

bool NacaAirfoil::contains_point(real_t x, real_t y, real_t lead_x, real_t lead_y) const {
    // Transform query point into local unrotated chord frame
    real_t dx = x - lead_x;
    real_t dy = y - lead_y;

    real_t alpha_rad = -aoa_deg_ * (PI / 180.0f);
    real_t cos_a = std::cos(-alpha_rad);
    real_t sin_a = std::sin(-alpha_rad);

    real_t x_chord = dx * cos_a - dy * sin_a;
    real_t y_chord = dx * sin_a + dy * cos_a;

    if (x_chord < 0.0f || x_chord > chord_) {
        return false;
    }

    real_t xc = x_chord / chord_;
    real_t yt = 5.0f * thickness_ * (
        0.2969f * std::sqrt(xc) -
        0.1260f * xc -
        0.3516f * xc * xc +
        0.2843f * xc * xc * xc -
        0.1036f * xc * xc * xc * xc
    ) * chord_;

    real_t yc = 0.0f;
    if (camber_pos_ > 1e-6f) {
        if (xc < camber_pos_) {
            yc = (camber_ / (camber_pos_ * camber_pos_)) * (2.0f * camber_pos_ * xc - xc * xc) * chord_;
        } else {
            real_t one_minus_p = 1.0f - camber_pos_;
            yc = (camber_ / (one_minus_p * one_minus_p)) * ((1.0f - 2.0f * camber_pos_) + 2.0f * camber_pos_ * xc - xc * xc) * chord_;
        }
    }

    return std::abs(y_chord - yc) <= yt;
}

void NacaAirfoil::rasterize(LatticeGrid2D& grid, real_t lead_x, real_t lead_y) const {
    real_t rad_bound = chord_ * 1.5f;
    int min_x = std::max(0, static_cast<int>(lead_x - rad_bound));
    int max_x = std::min(grid.nx() - 1, static_cast<int>(lead_x + rad_bound));
    int min_y = std::max(0, static_cast<int>(lead_y - rad_bound));
    int max_y = std::min(grid.ny() - 1, static_cast<int>(lead_y + rad_bound));

    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            if (contains_point(static_cast<real_t>(x), static_cast<real_t>(y), lead_x, lead_y)) {
                grid.set_cell_type(x, y, CellType::Solid);
            }
        }
    }
}

void NacaAirfoil::voxelize_extruded(
    LatticeGrid3D& grid,
    real_t lead_x,
    real_t lead_y,
    int z_start,
    int z_end
) const {
    real_t rad_bound = chord_ * 1.5f;
    int min_x = std::max(0, static_cast<int>(lead_x - rad_bound));
    int max_x = std::min(grid.nx() - 1, static_cast<int>(lead_x + rad_bound));
    int min_y = std::max(0, static_cast<int>(lead_y - rad_bound));
    int max_y = std::min(grid.ny() - 1, static_cast<int>(lead_y + rad_bound));

    int z0 = std::max(0, z_start);
    int z1 = std::min(grid.nz() - 1, z_end);

    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            if (contains_point(static_cast<real_t>(x), static_cast<real_t>(y), lead_x, lead_y)) {
                for (int z = z0; z <= z1; ++z) {
                    grid.set_cell_type(x, y, z, CellType::Solid);
                }
            }
        }
    }
}

} // namespace wake
