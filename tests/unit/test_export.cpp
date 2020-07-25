#include "wake/export/terminal_renderer.hpp"
#include "wake/export/vtk_exporter.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>

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

void test_vtk_exporter() {
    std::cout << "Testing VTK XML ImageData exporter..." << std::endl;

    LatticeGrid2D grid2d(32, 16, 0.02f);
    grid2d.initialize(1.0f, 0.04f, 0.01f);
    const std::string vti_2d_path = "test_output_2d.vti";
    bool ok2d = VtkExporter::export_vti_2d(vti_2d_path, grid2d);
    assert(ok2d);

    std::ifstream check2d(vti_2d_path);
    assert(check2d.is_open());
    std::string line;
    bool found_vti = false;
    while (std::getline(check2d, line)) {
        if (line.find("<VTKFile type=\"ImageData\"") != std::string::npos) {
            found_vti = true;
            break;
        }
    }
    assert(found_vti);
    check2d.close();
    std::remove(vti_2d_path.c_str());

    // 3D VTI export test
    LatticeGrid3D grid3d(16, 16, 16, 0.02f);
    grid3d.initialize(1.0f, 0.02f, 0.0f, 0.0f);
    const std::string vti_3d_path = "test_output_3d.vti";
    bool ok3d = VtkExporter::export_vti_3d(vti_3d_path, grid3d);
    assert(ok3d);
    std::remove(vti_3d_path.c_str());

    std::cout << "VTK XML exporter PASSED." << std::endl;
}

int main() {
    test_terminal_renderer();
    test_vtk_exporter();
    std::cout << "All export tests passed successfully!" << std::endl;
    return 0;
}
