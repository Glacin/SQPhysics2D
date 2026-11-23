cmake_minimum_required(VERSION 3.24)

project(sqphysics2d
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SQP_BUILD_TESTS "Build tests" ON)
option(SQP_BUILD_EXAMPLES "Build examples" ON)
option(SQP_DEBUG_VALIDATION "Runtime debug validations" ON)
option(SQP_ENABLE_PERFORMANCE_PROFILING "Enable runtime performance profiling timers" ON)

option(BUILD_SHARED_LIBS "Build shared library" OFF)

set(SQP_PRECISION "float" CACHE STRING "Default real_t precision")
set_property(CACHE SQP_PRECISION PROPERTY STRINGS float double)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

add_subdirectory(src)

if(SQP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(SQP_BUILD_TESTS AND EXISTS "${PROJECT_SOURCE_DIR}/tests")
    enable_testing()

    # Local vendored Catch2
    add_subdirectory(extern/catch2)

    add_subdirectory(tests)
endif()

install(
    TARGETS sqphysics2d
    EXPORT sqphysics2dTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(
    FILES
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/build_config.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/settings.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/material.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/shape.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/body.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/joint.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/world.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/events.h"
    DESTINATION include/sqphysics2d
)

install(
    FILES
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/types.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/math.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/contact.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/contact_solver.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/constraint_graph.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/data_storage.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/core_world.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/broad_phase.h"
        "${PROJECT_SOURCE_DIR}/include/sqphysics2d/detail/shape_query_view.h"
    DESTINATION include/sqphysics2d/detail
)

install(
    EXPORT sqphysics2dTargets
    NAMESPACE sqp::
    DESTINATION lib/cmake/sqphysics2d
)

install(
    FILES cmake/sqphysics2dConfig.cmake
    DESTINATION lib/cmake/sqphysics2d
)
