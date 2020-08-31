add_library(wake_runner STATIC
    runner.cpp
)

target_include_directories(wake_runner PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(wake_runner PUBLIC
    wake_core
    wake_scenarios
)
