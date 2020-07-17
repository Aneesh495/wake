#include "wake/geometry/primitives.hpp"
#include "wake/geometry/obj_loader.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace wake;

void test_primitives_2d() {
    std::cout << "Testing 2D geometry primitives..." << std::endl;

    LatticeGrid2D grid(64, 64, 0.02f);
    grid.initialize();

    GeometryPrimitives::rasterize_circle(grid, 32.0f, 32.0f, 10.0f);
    assert(grid.get_cell_type(32, 32) == CellType::Solid);
    assert(grid.get_cell_type(0, 0) == CellType::Fluid);

    GeometryPrimitives::rasterize_rectangle(grid, 5.0f, 5.0f, 15.0f, 15.0f);
    assert(grid.get_cell_type(10, 10) == CellType::Solid);

    GeometryPrimitives::rasterize_rotated_plate(grid, 45.0f, 45.0f, 20.0f, 4.0f, 45.0f);
    assert(grid.get_cell_type(45, 45) == CellType::Solid);

    std::cout << "2D geometry primitives PASSED." << std::endl;
}

void test_obj_loader_and_voxelizer() {
    std::cout << "Testing Wavefront OBJ parsing and 3D voxelization..." << std::endl;

    // Unit cube mesh OBJ string
    std::string cube_obj =
        "v -1.0 -1.0 -1.0\n"
        "v  1.0 -1.0 -1.0\n"
        "v  1.0  1.0 -1.0\n"
        "v -1.0  1.0 -1.0\n"
        "v -1.0 -1.0  1.0\n"
        "v  1.0 -1.0  1.0\n"
        "v  1.0  1.0  1.0\n"
        "v -1.0  1.0  1.0\n"
        "f 1 2 3 4\n"
        "f 5 8 7 6\n"
        "f 1 5 6 2\n"
        "f 2 6 7 3\n"
        "f 3 7 8 4\n"
        "f 5 1 4 8\n";

    ObjLoader loader;
    bool loaded = loader.load_from_string(cube_obj);
    assert(loaded);
    assert(loader.vertex_count() == 8);
    // 6 quad faces triangulate to 12 triangles
    assert(loader.triangle_count() == 12);

    LatticeGrid3D grid(32, 32, 32, 0.02f);
    grid.initialize();

    // Scale and fit cube into center of domain [10, 10, 10] to [22, 22, 22]
    loader.fit_into_box(Vec3(10.0f, 10.0f, 10.0f), Vec3(22.0f, 22.0f, 22.0f));
    loader.voxelize_into_grid(grid, true);

    // Center of cube must be solid
    assert(grid.get_cell_type(16, 16, 16) == CellType::Solid);

    // Far corners outside cube must remain fluid
    assert(grid.get_cell_type(2, 2, 2) == CellType::Fluid);
    assert(grid.get_cell_type(30, 30, 30) == CellType::Fluid);

    std::cout << "Wavefront OBJ loader and voxelizer PASSED." << std::endl;
}

int main() {
    test_primitives_2d();
    test_obj_loader_and_voxelizer();
    std::cout << "All geometry tests completed successfully!" << std::endl;
    return 0;
}
