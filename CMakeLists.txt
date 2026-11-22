cmake_minimum_required(VERSION 3.25)

project(litho_fdtd VERSION 0.1.0 DESCRIPTION "C++23 3D FDTD mask-spectrum solver" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LITHO_FDTD_BUILD_TESTS "Build litho_fdtd tests" ${PROJECT_IS_TOP_LEVEL})
option(LITHO_FDTD_BUILD_EXAMPLES "Build litho_fdtd examples" ${PROJECT_IS_TOP_LEVEL})
option(LITHO_FDTD_BUILD_BENCHMARKS "Build the throughput benchmark" OFF)
option(LITHO_FDTD_WARNINGS_AS_ERRORS "Treat project warnings as errors" OFF)

if(PROJECT_IS_TOP_LEVEL AND CMAKE_EXPORT_COMPILE_COMMANDS
   AND CMAKE_GENERATOR MATCHES "Ninja|Makefiles")
    set(_litho_compile_commands "${PROJECT_SOURCE_DIR}/build/compile_commands.json")
    if(NOT "${PROJECT_BINARY_DIR}/compile_commands.json" STREQUAL "${_litho_compile_commands}")
        file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/build")
        file(CREATE_LINK "${PROJECT_BINARY_DIR}/compile_commands.json"
            "${_litho_compile_commands}" SYMBOLIC)
    endif()
endif()

if(NOT TARGET litho_geometry::litho_geometry)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../litho_geometry"
        "${CMAKE_CURRENT_BINARY_DIR}/litho_geometry")
endif()

add_library(litho_fdtd
    src/mask.cpp
    src/structure.cpp
    src/spectrum.cpp
    src/solver.cpp
)
add_library(litho_fdtd::litho_fdtd ALIAS litho_fdtd)
target_compile_features(litho_fdtd PUBLIC cxx_std_23)
target_include_directories(litho_fdtd PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(litho_fdtd PUBLIC litho_geometry::litho_geometry)

if(NOT TARGET parallel)
    # The sibling parallel scheduler is the only worker engine. Disable its
    # network-fetched tests when vendored as a subproject.
    set(_litho_fdtd_saved_build_testing "${BUILD_TESTING}")
    set(BUILD_TESTING OFF)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../parallel"
        "${CMAKE_CURRENT_BINARY_DIR}/parallel" EXCLUDE_FROM_ALL)
    set(BUILD_TESTING "${_litho_fdtd_saved_build_testing}")
endif()
target_link_libraries(litho_fdtd PRIVATE parallel)

if(MSVC)
    target_compile_options(litho_fdtd PRIVATE /W4)
else()
    target_compile_options(litho_fdtd PRIVATE -Wall -Wextra -Wpedantic)
endif()
if(LITHO_FDTD_WARNINGS_AS_ERRORS)
    if(MSVC)
        target_compile_options(litho_fdtd PRIVATE /WX)
    else()
        target_compile_options(litho_fdtd PRIVATE -Werror)
    endif()
endif()

if(LITHO_FDTD_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(LITHO_FDTD_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(LITHO_FDTD_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
