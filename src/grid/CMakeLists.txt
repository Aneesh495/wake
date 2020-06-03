add_library(wake_grid STATIC
    lattice_grid_2d.cpp
    lattice_grid_3d.cpp
)

target_include_directories(wake_grid PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(wake_grid PUBLIC
    wake_core
    wake_ispc_kernels
)
