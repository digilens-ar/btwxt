/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

// Standard
#include <chrono>
#include <iostream>

// vendor
#include <gmock/gmock.h>
#include "gmock/gmock-matchers.h"
#include <gtest/gtest.h>

// btwxt
#include <btwxt/btwxt.h>
#include "fixtures/public-fixtures.h"
#include "spdlog/spdlog.h"

namespace Btwxt {

TEST_F(FunctionFixture, scipy_3d_grid)
{
    // Based on
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridinterpolator.value().html
    grid = {
        GridAxis(linspace(1, 4, 11)),
        GridAxis(linspace(4, 7, 22)),
        GridAxis(linspace(7, 9, 33))
    };

    functions = {[](std::vector<double> x) -> double {
        return 2 * x[0] * x[0] * x[0] + 3 * x[1] * x[1] - x[2];
    }};
    setup();

    const double epsilon = 0.0001;
    double result;
    double expected_value;

    target = {2.1, 6.2, 8.3};
    result = interpolator->solve(target)[0];
    expected_value = 125.80469388; // Interpolated value from example
    EXPECT_NEAR(result, expected_value, epsilon);

    target = {3.3, 5.2, 7.1};
    result = interpolator->solve(target)[0];
    expected_value = 146.30069388; // Interpolated value from example
    EXPECT_NEAR(result, expected_value, epsilon);
}


TEST_F(GridFixture, four_point_1d_cubic_interpolate)
{

    grid = {GridAxis({0, 2, 5, 10}, InterpolationMethod::cubic)};
    data_sets = {{6, 5, 4, 3}};
    target = {2.5};
    setup();

    const double expected_value = 4.804398;
    const double epsilon = 0.0001;
    double result = interpolator->solve(target)[0];
    EXPECT_NEAR(result, expected_value, epsilon);
}

TEST_F(GridFixture, empty_grid_throw_test)
{
    grid = {{}};
    data_sets = {{}};
    target = {};
    EXPECT_THROW(setup(), std::runtime_error);
}

TEST_F(GridFixture, two_point_cubic_1d_interpolate)
{
    grid = {GridAxis({0, 10}, InterpolationMethod::cubic)};
    data_sets = {{6, 3}};
    target = {2.5};
    setup();
    double result = interpolator->solve(target)[0];

    EXPECT_NEAR(result, 5.25, 0.0001);
}

TEST_F(GridFixture, get_neighboring_indices)
{
    grid = {GridAxis({0, 1, 2}), GridAxis({0, 1, 2})};

    // data_set[i] = i useful for testing
    // clang-format off
    data_sets = {{
    //  0  1  2 < dim 2
        0, 1, 2, // 0 dim 1
        3, 4, 5, // 1  "
        6, 7, 8  // 2  "
    }};
    // clang-format on
    setup();
}


TEST_F(Grid2DFixture, interpolate)
{
    auto result = interpolator.value().solve(target);

    // All values, current target
    EXPECT_THAT(result, testing::ElementsAre(testing::DoubleEq(4.2), testing::DoubleEq(8.4)));
    // Single value, current target
    double d_result = result[0];
    EXPECT_DOUBLE_EQ(d_result, 4.2);

    std::vector<double> another_target = {8.1, 4.2};
    // All values, fresh target
    result = interpolator.value().solve(another_target);
    EXPECT_THAT(result, testing::ElementsAre(testing::DoubleEq(3.189), testing::DoubleEq(6.378)));
    // Single value, fresh target
    d_result = result[1];
    EXPECT_DOUBLE_EQ(d_result, 6.378);
}

TEST_F(Grid2DFixture, extrapolate)
{
    // axis1 is designated constant extrapolation
    target = {10, 3};
    EXPECT_THROW(interpolator->solve(target), std::runtime_error);
}

TEST_F(Grid2DFixtureCubic, cubic_interpolate)
{
    auto result = interpolator.value().solve(target);
    // All values, current target
    EXPECT_THAT(result, testing::ElementsAre(testing::DoubleEq(4.416), testing::DoubleEq(8.832)));
}

TEST_F(Function4DFixture, calculate)
{
    auto result = interpolator.value().solve(target);
    EXPECT_NEAR(result[0], functions[0](target), 0.02);
    EXPECT_DOUBLE_EQ(result[1], functions[1](target));
}

TEST_F(Function4DFixture, timer)
{
    auto result = interpolator.value().solve(target);

    // Get starting time point
    auto start = std::chrono::high_resolution_clock::now();
    // Get ending time point
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    spdlog::info("Time taken by interpolation: {} microseconds", duration.count());

    // time running the functions straight
    start = std::chrono::high_resolution_clock::now();
    // double r0 = fn0(target[0], target[1], target[2], target[3]);
    // double r1 = fn1(target[0], target[1], target[2], target[3]);
    // Get ending time point
    stop = std::chrono::high_resolution_clock::now();
    auto nano_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
   spdlog::info(
        fmt::format("Time taken by direct functions: {} nanoseconds", nano_duration.count()));
}

TEST_F(Function4DFixture, multi_timer)
{
    std::vector<std::vector<double>> set_of_targets = {{0.1, 0.1, 0.1, 0.1},
                                                       {3.3, 2.2, 4.1, 1.4},
                                                       {2.1, 1.6, 1.6, 2.1},
                                                       {3.7, 4.3, 0.8, 2.1},
                                                       {1.9, 3.4, 1.2, 1.1},
                                                       {3.3, 3.8, 1.6, 3.0},
                                                       {0.3, 1.0, 2.4, 1.1},
                                                       {3.1, 1.9, 2.9, 3.3},
                                                       {4.2, 2.7, 1.3, 4.4},
                                                       {2.1, 2.9, 1.8, 1.9}};

    for (std::size_t count = 0; count < 10; count++) {
        // Get starting time point
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto& target : set_of_targets) {
            auto result = interpolator.value().solve(target);
        }
        // Get ending time point
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        spdlog::info(
            "Time taken by ten interpolations: {} microseconds", duration.count());
    }
}

TEST(GridPointDataSet, wrong_size)
{
    EXPECT_THROW(RegularGridInterpolator({GridAxis({1.})}, {GridPointDataSet({1., 1.})}),
                     std::runtime_error);
}

TEST(CartesianProduct, cartesian_product)
{
    std::vector<std::vector<short>> result =
        cartesian_product(std::vector<std::vector<short>> {{1, 2, 3}, {4, 5}, {6, 7, 8, 9}});
    EXPECT_EQ(result.size(), 3u * 2u * 4u);
    EXPECT_THAT(result[0], testing::ElementsAre(1, 4, 6));
    EXPECT_THAT(result[1], testing::ElementsAre(1, 4, 7));
    EXPECT_THAT(result[10], testing::ElementsAre(2, 4, 8));
    EXPECT_THAT(result[3 * 2 * 4 - 1], testing::ElementsAre(3, 5, 9));
}

} // namespace Btwxt