add_library(wake_scenarios STATIC
    scenario_registry.cpp
    karman_vortex.cpp
    naca_airfoil_stall.cpp
    lid_driven_cavity.cpp
    sphere_3d_wake.cpp
)

target_include_directories(wake_scenarios PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(wake_scenarios PUBLIC
    wake_core
    wake_grid
    wake_forces
    wake_geometry
    wake_export
    wake_analytics
)
