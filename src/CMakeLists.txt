
add_library(btwxt
    regular-grid-interpolator-implementation.h
    regular-grid-interpolator-implementation.cpp
    regular-grid-interpolator.cpp
    grid-axis.cpp
)

target_include_directories(btwxt PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(btwxt PRIVATE spdlog::spdlog)