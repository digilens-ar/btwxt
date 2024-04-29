
add_library(btwxt
    STATIC
    regular-grid-interpolator-implementation.h
    regular-grid-interpolator-implementation.cpp
    regular-grid-interpolator.cpp
    grid-axis.cpp
)

target_link_libraries(btwxt PRIVATE spdlog::spdlog)
target_sources(btwxt PUBLIC
    FILE_SET public_headers
    TYPE HEADERS
    BASE_DIRS "${PROJECT_SOURCE_DIR}/include"
    FILES "${PROJECT_SOURCE_DIR}/include/btwxt/btwxt.h" "${PROJECT_SOURCE_DIR}/include/btwxt/grid-axis.h" "${PROJECT_SOURCE_DIR}/include/btwxt/grid-point-data.h" "${PROJECT_SOURCE_DIR}/include/btwxt/regular-grid-interpolator.h"
)
install(TARGETS btwxt FILE_SET public_headers)