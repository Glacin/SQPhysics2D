file(GLOB SQP_BENCHMARK_SOURCES CONFIGURE_DEPENDS "benchmark/*.cpp")

add_executable(sqphysics2d_benchmark ${SQP_BENCHMARK_SOURCES})
target_link_libraries(sqphysics2d_benchmark PRIVATE sqphysics2d)
target_include_directories(sqphysics2d_benchmark PRIVATE "${PROJECT_SOURCE_DIR}/src")
set_target_properties(sqphysics2d_benchmark PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

add_executable(sqphysics2d_replay replay.cpp)
target_link_libraries(sqphysics2d_replay PRIVATE sqphysics2d)
set_target_properties(sqphysics2d_replay PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
