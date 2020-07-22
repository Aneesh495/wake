#pragma once

#include "wake/core/types.hpp"
#include "wake/grid/lattice_grid_2d.hpp"
#include "wake/grid/lattice_grid_3d.hpp"
#include <string>
#include <vector>
#include <ostream>

namespace wake {

enum class Colormap {
    Viridis,
    Inferno,
    CoolWarm,
    Jet
};

struct RGBColor {
    uint8_t r{0};
    uint8_t g{0};
    uint8_t b{0};
};

class TerminalRenderer {
public:
    explicit TerminalRenderer(
        int term_width = 80,
        int term_height = 24,
        Colormap cmap = Colormap::Viridis
    );

    void set_colormap(Colormap cmap) noexcept { colormap_ = cmap; }
    void set_dimensions(int width, int height) noexcept { width_ = width; height_ = height; }

    // Map normalized value [0, 1] to RGBColor
    RGBColor map_value(real_t t) const;

    // Render 2D scalar field (or speed |u|) using half-block characters
    std::string render_2d(
        const LatticeGrid2D& grid,
        const real_t* scalar_field = nullptr, // if null, uses velocity magnitude
        real_t min_val = 0.0f,
        real_t max_val = 0.15f
    ) const;

    // Render 2D vorticity field with diverging CoolWarm colormap centered at zero
    std::string render_vorticity_2d(
        const LatticeGrid2D& grid,
        const real_t* vorticity,
        real_t max_vort = 0.05f
    ) const;

    // Render 2D mid-plane slice of 3D grid
    std::string render_slice_3d(
        const LatticeGrid3D& grid,
        int z_slice = -1, // if -1, mid-plane nz / 2
        real_t min_val = 0.0f,
        real_t max_val = 0.15f
    ) const;

    // Print frame to standard output with cursor home
    void display(std::ostream& out, const std::string& frame) const;

private:
    int width_{80};
    int height_{24};
    Colormap colormap_{Colormap::Viridis};

    static RGBColor sample_colormap_viridis(real_t t);
    static RGBColor sample_colormap_inferno(real_t t);
    static RGBColor sample_colormap_coolwarm(real_t t);
    static RGBColor sample_colormap_jet(real_t t);
};

} // namespace wake
