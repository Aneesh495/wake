#include "wake/export/terminal_renderer.hpp"
#include <cmath>
#include <algorithm>
#include <sstream>

namespace wake {

TerminalRenderer::TerminalRenderer(int term_width, int term_height, Colormap cmap)
    : width_(term_width), height_(term_height), colormap_(cmap) {}

RGBColor TerminalRenderer::sample_colormap_viridis(real_t t) {
    t = std::clamp(t, 0.0f, 1.0f);
    // Piecewise cubic polynomial approximation to Matplotlib Viridis
    real_t r = -0.019f + t * (1.78f + t * (-3.29f + t * 2.36f));
    real_t g =  0.005f + t * (0.83f + t * ( 0.94f - t * 0.88f));
    real_t b =  0.284f + t * (1.08f + t * (-2.61f + t * 1.55f));

    r = std::clamp(r, 0.0f, 1.0f);
    g = std::clamp(g, 0.0f, 1.0f);
    b = std::clamp(b, 0.0f, 1.0f);

    return {
        static_cast<uint8_t>(r * 255.0f),
        static_cast<uint8_t>(g * 255.0f),
        static_cast<uint8_t>(b * 255.0f)
    };
}

RGBColor TerminalRenderer::sample_colormap_inferno(real_t t) {
    t = std::clamp(t, 0.0f, 1.0f);
    real_t r = 0.001f + t * (2.81f + t * (-3.42f + t * 1.61f));
    real_t g = 0.001f + t * (0.05f + t * ( 1.25f - t * 0.30f));
    real_t b = 0.015f + t * (1.43f - t * ( 2.80f - t * 1.35f));

    r = std::clamp(r, 0.0f, 1.0f);
    g = std::clamp(g, 0.0f, 1.0f);
    b = std::clamp(b, 0.0f, 1.0f);

    return {
        static_cast<uint8_t>(r * 255.0f),
        static_cast<uint8_t>(g * 255.0f),
        static_cast<uint8_t>(b * 255.0f)
    };
}

RGBColor TerminalRenderer::sample_colormap_coolwarm(real_t t) {
    // Diverging blue -> white -> red
    t = std::clamp(t, 0.0f, 1.0f);
    real_t r, g, b;
    if (t < 0.5f) {
        real_t s = t * 2.0f; // 0 to 1
        r = 0.23f + 0.77f * s;
        g = 0.29f + 0.71f * s;
        b = 0.75f + 0.25f * s;
    } else {
        real_t s = (t - 0.5f) * 2.0f; // 0 to 1
        r = 1.00f - 0.30f * s;
        g = 1.00f - 0.75f * s;
        b = 1.00f - 0.75f * s;
    }
    return {
        static_cast<uint8_t>(std::clamp(r, 0.0f, 1.0f) * 255.0f),
        static_cast<uint8_t>(std::clamp(g, 0.0f, 1.0f) * 255.0f),
        static_cast<uint8_t>(std::clamp(b, 0.0f, 1.0f) * 255.0f)
    };
}

RGBColor TerminalRenderer::sample_colormap_jet(real_t t) {
    t = std::clamp(t, 0.0f, 1.0f);
    real_t r = std::clamp(1.5f - std::abs(4.0f * t - 3.0f), 0.0f, 1.0f);
    real_t g = std::clamp(1.5f - std::abs(4.0f * t - 2.0f), 0.0f, 1.0f);
    real_t b = std::clamp(1.5f - std::abs(4.0f * t - 1.0f), 0.0f, 1.0f);
    return {
        static_cast<uint8_t>(r * 255.0f),
        static_cast<uint8_t>(g * 255.0f),
        static_cast<uint8_t>(b * 255.0f)
    };
}

RGBColor TerminalRenderer::map_value(real_t t) const {
    switch (colormap_) {
        case Colormap::Inferno: return sample_colormap_inferno(t);
        case Colormap::CoolWarm: return sample_colormap_coolwarm(t);
        case Colormap::Jet: return sample_colormap_jet(t);
        case Colormap::Viridis:
        default: return sample_colormap_viridis(t);
    }
}

std::string TerminalRenderer::render_2d(
    const LatticeGrid2D& grid,
    const real_t* scalar_field,
    real_t min_val,
    real_t max_val
) const {
    std::ostringstream ss;
    int nx = grid.nx();
    int ny = grid.ny();

    real_t val_range = (max_val - min_val > 1e-8f) ? (max_val - min_val) : 1.0f;
    int pixel_rows = height_ * 2;

    const RGBColor solid_color{110, 110, 115};

    auto sample_point = [&](int px, int py) -> RGBColor {
        int gx = std::clamp(px * nx / width_, 0, nx - 1);
        int gy = std::clamp(py * ny / pixel_rows, 0, ny - 1);

        if (grid.get_cell_type(gx, gy) == CellType::Solid) {
            return solid_color;
        }

        real_t val = 0.0f;
        if (scalar_field) {
            val = scalar_field[gy * nx + gx];
        } else {
            val = grid.speed(gx, gy);
        }

        real_t norm = (val - min_val) / val_range;
        return map_value(norm);
    };

    // Render using Unicode upper half block: U+2580 ('▀')
    // Top pixel = foreground, Bottom pixel = background
    for (int row = height_ - 1; row >= 0; --row) {
        int py_top = row * 2 + 1;
        int py_bot = row * 2;

        for (int col = 0; col < width_; ++col) {
            RGBColor top_c = sample_point(col, py_top);
            RGBColor bot_c = sample_point(col, py_bot);

            ss << "\033[38;2;" << static_cast<int>(top_c.r) << ";"
               << static_cast<int>(top_c.g) << ";"
               << static_cast<int>(top_c.b) << "m"
               << "\033[48;2;" << static_cast<int>(bot_c.r) << ";"
               << static_cast<int>(bot_c.g) << ";"
               << static_cast<int>(bot_c.b) << "m"
               << "\u2580";
        }
        ss << "\033[0m\n";
    }

    return ss.str();
}

std::string TerminalRenderer::render_vorticity_2d(
    const LatticeGrid2D& grid,
    const real_t* vorticity,
    real_t max_vort
) const {
    std::ostringstream ss;
    int nx = grid.nx();
    int ny = grid.ny();

    int pixel_rows = height_ * 2;
    real_t denom = (max_vort > 1e-8f) ? (2.0f * max_vort) : 1.0f;
    const RGBColor solid_color{110, 110, 115};

    auto sample_vort = [&](int px, int py) -> RGBColor {
        int gx = std::clamp(px * nx / width_, 0, nx - 1);
        int gy = std::clamp(py * ny / pixel_rows, 0, ny - 1);

        if (grid.get_cell_type(gx, gy) == CellType::Solid) {
            return solid_color;
        }

        real_t omega = vorticity[gy * nx + gx];
        real_t norm = (omega + max_vort) / denom;
        return sample_colormap_coolwarm(norm);
    };

    for (int row = height_ - 1; row >= 0; --row) {
        int py_top = row * 2 + 1;
        int py_bot = row * 2;

        for (int col = 0; col < width_; ++col) {
            RGBColor top_c = sample_vort(col, py_top);
            RGBColor bot_c = sample_vort(col, py_bot);

            ss << "\033[38;2;" << static_cast<int>(top_c.r) << ";"
               << static_cast<int>(top_c.g) << ";"
               << static_cast<int>(top_c.b) << "m"
               << "\033[48;2;" << static_cast<int>(bot_c.r) << ";"
               << static_cast<int>(bot_c.g) << ";"
               << static_cast<int>(bot_c.b) << "m"
               << "\u2580";
        }
        ss << "\033[0m\n";
    }

    return ss.str();
}

std::string TerminalRenderer::render_slice_3d(
    const LatticeGrid3D& grid,
    int z_slice,
    real_t min_val,
    real_t max_val
) const {
    std::ostringstream ss;
    int nx = grid.nx();
    int ny = grid.ny();
    int nz = grid.nz();

    if (z_slice < 0 || z_slice >= nz) {
        z_slice = nz / 2;
    }

    real_t val_range = (max_val - min_val > 1e-8f) ? (max_val - min_val) : 1.0f;
    int pixel_rows = height_ * 2;
    const RGBColor solid_color{110, 110, 115};

    auto sample_point = [&](int px, int py) -> RGBColor {
        int gx = std::clamp(px * nx / width_, 0, nx - 1);
        int gy = std::clamp(py * ny / pixel_rows, 0, ny - 1);

        if (grid.get_cell_type(gx, gy, z_slice) == CellType::Solid) {
            return solid_color;
        }

        real_t val = grid.speed(gx, gy, z_slice);
        real_t norm = (val - min_val) / val_range;
        return map_value(norm);
    };

    for (int row = height_ - 1; row >= 0; --row) {
        int py_top = row * 2 + 1;
        int py_bot = row * 2;

        for (int col = 0; col < width_; ++col) {
            RGBColor top_c = sample_point(col, py_top);
            RGBColor bot_c = sample_point(col, py_bot);

            ss << "\033[38;2;" << static_cast<int>(top_c.r) << ";"
               << static_cast<int>(top_c.g) << ";"
               << static_cast<int>(top_c.b) << "m"
               << "\033[48;2;" << static_cast<int>(bot_c.r) << ";"
               << static_cast<int>(bot_c.g) << ";"
               << static_cast<int>(bot_c.b) << "m"
               << "\u2580";
        }
        ss << "\033[0m\n";
    }

    return ss.str();
}

void TerminalRenderer::display(std::ostream& out, const std::string& frame) const {
    out << "\033[H" << frame << std::flush;
}

} // namespace wake
