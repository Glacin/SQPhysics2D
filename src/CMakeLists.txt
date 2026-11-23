file(GLOB_RECURSE SQP_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_library(sqphysics2d
    ${SQP_SOURCES}
)

target_include_directories(sqphysics2d
    PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:include>"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_features(sqphysics2d
    PUBLIC
        cxx_std_20
)

target_compile_definitions(sqphysics2d
    PRIVATE
        SQP_BUILD
)

if (SQP_DEBUG_VALIDATION)
    target_compile_definitions(sqphysics2d
        PUBLIC 
            SQP_DEBUG_VALIDATION
)
endif()

if(SQP_ENABLE_PERFORMANCE_PROFILING)
    target_compile_definitions(sqphysics2d PUBLIC SQP_ENABLE_PERFORMANCE_PROFILING)
endif()

if(SQP_PRECISION STREQUAL "double")
    target_compile_definitions(sqphysics2d
        PUBLIC
            SQP_REAL_DOUBLE
    )
elseif(SQP_PRECISION STREQUAL "float")
    target_compile_definitions(sqphysics2d
        PUBLIC
            SQP_REAL_FLOAT
    )
else()
    message(FATAL_ERROR "Invalid SQP_PRECISION: ${SQP_PRECISION}. Expected float or double.")
endif()

if(MSVC)
    target_compile_options(sqphysics2d
        PRIVATE
            /W4
            /permissive-
    )
else()
    target_compile_options(sqphysics2d
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
    )
endif()

target_compile_definitions(sqphysics2d
    PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
)
