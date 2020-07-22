#include "wake/export/terminal_renderer.hpp"
#include <iostream>
#include <cassert>

using namespace wake;

void test_terminal_renderer() {
    std::cout << "Testing ANSI TrueColor terminal renderer..." << std::endl;

    LatticeGrid2D grid2d(64, 32, 0.02f);
    grid2d.initialize(1.0f, 0.05f, 0.0f);
    grid2d.set_solid_circle(16.0f, 16.0f, 5.0f);

    TerminalRenderer renderer(40, 16, Colormap::Viridis);
    std::string frame = renderer.render_2d(grid2d);

    assert(!frame.empty());
    // Must contain ANSI escape sequence \033[
    assert(frame.find("\033[") != std::string::npos);

    // Test 3D slice rendering
    LatticeGrid3D grid3d(32, 32, 16, 0.02f);
    grid3d.initialize(1.0f, 0.05f, 0.0f, 0.0f);
    std::string frame3d = renderer.render_slice_3d(grid3d, 8);
    assert(!frame3d.empty());

    std::cout << "Terminal renderer PASSED." << std::endl;
}

int main() {
    test_terminal_renderer();
    std::cout << "All export tests passed successfully!" << std::endl;
    return 0;
}
