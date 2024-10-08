
add_library(btwxt
    STATIC
    "regular-grid-interpolator.cpp"
    "grid-axis.cpp"
)

target_link_libraries(btwxt PRIVATE spdlog::spdlog)
target_sources(btwxt PUBLIC
    FILE_SET public_headers
    TYPE HEADERS
    BASE_DIRS "${PROJECT_SOURCE_DIR}/include"
    FILES 
    "${PROJECT_SOURCE_DIR}/include/btwxt/btwxt.hpp"
    "${PROJECT_SOURCE_DIR}/include/btwxt/grid-axis.hpp"
    "${PROJECT_SOURCE_DIR}/include/btwxt/regular-grid-interpolator.hpp"
    "${PROJECT_SOURCE_DIR}/include/btwxt/utility-p.hpp"
)

add_library(btwxt::btwxt ALIAS btwxt)  # Make this library accessible through the same name if you import or if you add subdirectory

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(
    TARGETS btwxt 
    EXPORT btwxt-export
    FILE_SET public_headers
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)


write_basic_package_version_file(
    btwxt-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(
    EXPORT btwxt-export
    NAMESPACE btwxt::
    FILE btwxt-targets.cmake
)

install(
    EXPORT btwxt-export
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/btwxt
    FILE btwxt-targets.cmake
    NAMESPACE btwxt::
)

install(
    FILES
        ${PROJECT_SOURCE_DIR}/cmake/btwxt-config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/btwxt-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/btwxt
)