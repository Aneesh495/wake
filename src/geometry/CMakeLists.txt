add_library(wake_geometry STATIC
    primitives.cpp
    obj_loader.cpp
)

target_include_directories(wake_geometry PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(wake_geometry PUBLIC
    wake_core
    wake_grid
)
