add_library(wake_ibm STATIC
    immersed_boundary.cpp
)

target_include_directories(wake_ibm PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(wake_ibm PUBLIC
    wake_core
    wake_ispc_kernels
)
