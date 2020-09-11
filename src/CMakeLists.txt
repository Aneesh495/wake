# Subsystems
add_subdirectory(core)
add_subdirectory(ispc)
add_subdirectory(grid)
add_subdirectory(forces)
add_subdirectory(ibm)
add_subdirectory(analytics)
add_subdirectory(geometry)
add_subdirectory(export)
add_subdirectory(thermal)
add_subdirectory(scenarios)
add_subdirectory(runner)

add_executable(wake main.cpp)
target_link_libraries(wake PRIVATE
    wake_runner
    wake_scenarios
    wake_thermal
    wake_export
    wake_geometry
    wake_analytics
    wake_forces
    wake_grid
    wake_core
    wake_ispc_kernels
)
