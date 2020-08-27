add_library(wake_thermal STATIC
    thermal_solver.cpp
)

target_include_directories(wake_thermal PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(wake_thermal PUBLIC
    wake_core
    wake_grid
    wake_ispc_kernels
)
