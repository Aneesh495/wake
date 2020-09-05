# Benchmarks CMakeLists.txt

add_executable(bench_d2q9_ispc_vs_scalar bench_d2q9_ispc_vs_scalar.cpp)
target_link_libraries(bench_d2q9_ispc_vs_scalar PRIVATE
    wake_core
    wake_ispc_kernels
)

add_executable(bench_d3q19_ispc_vs_scalar bench_d3q19_ispc_vs_scalar.cpp)
target_link_libraries(bench_d3q19_ispc_vs_scalar PRIVATE
    wake_core
    wake_ispc_kernels
)
