add_library(wake_export STATIC
    terminal_renderer.cpp
    vtk_exporter.cpp
    image_exporter.cpp
)

target_include_directories(wake_export PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(wake_export PUBLIC
    wake_core
    wake_grid
    wake_analytics
)
