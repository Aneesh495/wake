# Tests CMakeLists.txt

macro(wake_add_test test_name source_file)
    add_executable(${test_name} ${source_file})
    target_link_libraries(${test_name} PRIVATE
        wake_grid
        wake_core
        wake_ispc_kernels
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endmacro()

wake_add_test(test_d2q9_moments unit/test_d2q9_moments.cpp)
wake_add_test(test_d2q9_collision unit/test_d2q9_collision.cpp)
