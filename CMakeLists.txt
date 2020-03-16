cmake_minimum_required(VERSION 3.20)
project(wake LANGUAGES CXX ISPC)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Determine host architecture and configure ISPC compilation target
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    set(ISPC_TARGET "neon-i32x4")
    set(ISPC_ARCH "aarch64")
else()
    set(ISPC_TARGET "avx2-i32x8")
    set(ISPC_ARCH "x86-64")
endif()

message(STATUS "Configuring wake for architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Selected ISPC target: ${ISPC_TARGET}")

# Compiler flags
set(CMAKE_ISPC_FLAGS "${CMAKE_ISPC_FLAGS} -O3 --opt=fast-math --pic --target=${ISPC_TARGET}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pthread")

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Common directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/wake/ispc)

# Options
option(WAKE_BUILD_TESTS "Build wake unit, validation, and integration tests" ON)
option(WAKE_BUILD_BENCHMARKS "Build wake performance benchmarks" ON)

# Source lists will be populated modularly
add_subdirectory(src)

if(WAKE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(WAKE_BUILD_BENCHMARKS)
    add_subdirectory(bench)
endif()
