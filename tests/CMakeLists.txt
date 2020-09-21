# Tests CMakeLists.txt

macro(wake_add_test test_name source_file)
    add_executable(${test_name} ${source_file})
    target_link_libraries(${test_name} PRIVATE
        wake_runner
        wake_scenarios
        wake_thermal
        wake_export
        wake_geometry
        wake_analytics
        wake_forces
        wake_grid
        wake_core
        wake_ispc_kernels
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endmacro()

wake_add_test(test_d2q9_moments unit/test_d2q9_moments.cpp)
wake_add_test(test_d2q9_collision unit/test_d2q9_collision.cpp)
wake_add_test(test_poiseuille_convergence validation/test_poiseuille_convergence.cpp)
wake_add_test(test_d3q19_moments unit/test_d3q19_moments.cpp)
wake_add_test(test_d3q19_collision unit/test_d3q19_collision.cpp)
wake_add_test(test_taylor_green_decay validation/test_taylor_green_decay.cpp)
wake_add_test(test_forces_momentum_exchange unit/test_forces_momentum_exchange.cpp)
wake_add_test(test_field_analytics unit/test_field_analytics.cpp)
wake_add_test(test_geometry unit/test_geometry.cpp)
wake_add_test(test_export unit/test_export.cpp)
wake_add_test(test_scenarios unit/test_scenarios.cpp)
wake_add_test(test_runner unit/test_runner.cpp)
wake_add_test(test_karman_vortex_shedding integration/test_karman_vortex_shedding.cpp)
wake_add_test(test_full_pipeline integration/test_full_pipeline.cpp)
