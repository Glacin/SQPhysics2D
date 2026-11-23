cmake_minimum_required(VERSION 3.22)
project(sqp_box2d_pyramid_reference LANGUAGES C CXX)

set(BOX2D_SAMPLES OFF CACHE BOOL "" FORCE)
set(BOX2D_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(BOX2D_BENCHMARKS OFF CACHE BOOL "" FORCE)
add_subdirectory(../../../box2d box2d)

add_executable(sqp_box2d_pyramid_reference main.cpp)
target_link_libraries(sqp_box2d_pyramid_reference PRIVATE box2d::box2d)
