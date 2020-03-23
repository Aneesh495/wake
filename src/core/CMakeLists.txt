add_library(wake_core STATIC
    memory.cpp
)

target_include_directories(wake_core PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)
