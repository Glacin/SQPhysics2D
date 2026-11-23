cmake_minimum_required(VERSION 3.24)
project(sqphysics2d_installed_consumer LANGUAGES CXX)

find_package(sqphysics2d CONFIG REQUIRED)

add_executable(sqphysics2d_installed_replay ../../examples/replay.cpp)
target_link_libraries(sqphysics2d_installed_replay PRIVATE sqp::sqphysics2d)

enable_testing()
add_test(NAME installed_replay COMMAND sqphysics2d_installed_replay)
