#pragma once

#include "wake/core/types.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include <string>
#include <vector>

namespace wake {

struct AirfoilPoint {
    real_t x{0.0f};
    real_t y{0.0f};
};

class NacaAirfoil {
public:
    // Construct from 4-digit designation string, e.g. "0012", "2412", "4415"
    explicit NacaAirfoil(
        const std::string& code = "0012",
        real_t chord = 100.0f,
        real_t angle_of_attack_deg = 0.0f,
        size_t num_points = 200
    );

    // Reconfigure airfoil parameters
    void configure(
        const std::string& code,
        real_t chord,
        real_t angle_of_attack_deg,
        size_t num_points = 200
    );

    // Get polygon perimeter vertices in clockwise order (LE -> Upper -> TE -> Lower -> LE)
    const std::vector<AirfoilPoint>& polygon() const { return polygon_; }

    // Check if a point is inside the airfoil contour (in local or global frame)
    bool contains_point(real_t x, real_t y, real_t lead_x, real_t lead_y) const;

    // Rasterize 2D airfoil into LatticeGrid2D
    void rasterize(LatticeGrid2D& grid, real_t lead_x, real_t lead_y) const;

    // Extrude airfoil spanwise into LatticeGrid3D between z_start and z_end
    void voxelize_extruded(
        LatticeGrid3D& grid,
        real_t lead_x,
        real_t lead_y,
        int z_start,
        int z_end
    ) const;

    // Accessors
    real_t chord() const noexcept { return chord_; }
    real_t angle_of_attack() const noexcept { return aoa_deg_; }
    real_t max_thickness() const noexcept { return thickness_; }
    real_t max_camber() const noexcept { return camber_; }

private:
    void generate_coordinates(size_t num_points);

    std::string code_;
    real_t chord_{100.0f};
    real_t aoa_deg_{0.0f};
    real_t camber_{0.0f};      // m
    real_t camber_pos_{0.0f};  // p
    real_t thickness_{0.12f};  // t
    std::vector<AirfoilPoint> polygon_;
};

} // namespace wake
