
cmake_minimum_required(VERSION 3.10)

project(btwxt VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 23)

find_package(Git QUIET)
find_package(spdlog CONFIG REQUIRED)

option(BTWXT_BUILD_TESTING "Build BTWXT testing targets" OFF)

add_subdirectory(src)

if (BTWXT_BUILD_TESTING)
    find_package(GTest CONFIG REQUIRED)
    enable_testing()
    add_subdirectory(test)
endif ()
