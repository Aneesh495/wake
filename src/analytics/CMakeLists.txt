add_library(wake_analytics STATIC
    field_analytics.cpp
)

target_include_directories(wake_analytics PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(wake_analytics PUBLIC
    wake_grid
    wake_core
    wake_ispc_kernels
)
