# Test Executable
add_executable(btwxt_tests
    test-main.cpp
    fixtures/public-fixtures.h
    fixtures/implementation-fixtures.h
    btwxt-tests.cpp
    implementation-tests.cpp
    grid-axis-tests.cpp)

target_include_directories(btwxt_tests PRIVATE "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(btwxt_tests PRIVATE btwxt GTest::gtest GTest::gmock spdlog::spdlog)

include(GoogleTest)

gtest_discover_tests(${PROJECT_NAME}_tests TEST_PREFIX ${PROJECT_NAME}:)
