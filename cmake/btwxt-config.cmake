include(CMakeFindDependencyMacro)
find_dependency(spdlog)

include("${CMAKE_CURRENT_LIST_DIR}/btwxt-targets.cmake")

check_required_components(btwxt)