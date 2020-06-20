add_library(wake_forces STATIC
    aerodynamic_telemetry.cpp
    strouhal_analyzer.cpp
)

target_include_directories(wake_forces PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(wake_forces PUBLIC
    wake_core
    wake_ispc_kernels
)
