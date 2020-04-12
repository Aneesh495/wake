# ISPC Kernels Library
set(ISPC_SOURCES
    d2q9_moments.ispc
    d2q9_collision_bgk.ispc
    d2q9_collision_trt.ispc
)

add_library(wake_ispc_kernels STATIC
    ${ISPC_SOURCES}
)

target_include_directories(wake_ispc_kernels PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

set_target_properties(wake_ispc_kernels PROPERTIES
    ISPC_INSTRUCTION_SETS "${ISPC_TARGET}"
)
