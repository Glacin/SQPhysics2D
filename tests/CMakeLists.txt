file(GLOB_RECURSE SQP_TEST_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_executable(sqphysics2d_tests ${SQP_TEST_SOURCES})

target_include_directories(sqphysics2d_tests
    PRIVATE
        "${PROJECT_SOURCE_DIR}/src"
        "${PROJECT_SOURCE_DIR}/extern/catch2/extras"
)

target_link_libraries(sqphysics2d_tests
    PRIVATE
        sqphysics2d
        Catch2::Catch2WithMain
)

set_target_properties(sqphysics2d_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/extern/catch2/extras")
include(Catch)
catch_discover_tests(sqphysics2d_tests)

if(SQP_BUILD_EXAMPLES)
    add_test(NAME replay_example COMMAND sqphysics2d_replay)
    add_test(NAME benchmark_circle_box_bounce_ccd
             COMMAND sqphysics2d_benchmark circle_box_bounce --warmup 7 --ticks 1)
    add_test(NAME benchmark_determinism
             COMMAND sqphysics2d_benchmark circle_box_bounce --bodies 12
                     --warmup 2 --ticks 3 --check-determinism)
endif()
