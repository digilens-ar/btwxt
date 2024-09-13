/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

// Standard
#include <chrono>

// vendor
#include <fmt/format.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// btwxt
#include "fixtures/implementation-fixtures.h"

namespace Btwxt {

TEST_F(CubicImplementationFixture, spacing_multiplier)
{
    static constexpr std::size_t floor = 0;
    static constexpr std::size_t ceiling = 1;
    double result;
    result = interpolator.get_grid_axis(0).get_cubic_spacing_ratios(floor)[0];
    EXPECT_DOUBLE_EQ(result, 1.0);

    result = interpolator.get_grid_axis(0).get_cubic_spacing_ratios(ceiling)[0];
    EXPECT_DOUBLE_EQ(result, (10. - 6.) / (15.0 - 6.0));

    result = interpolator.get_grid_axis(0).get_cubic_spacing_ratios(floor)[1];
    EXPECT_DOUBLE_EQ(result, (15. - 10.) / (15.0 - 6.0));

    result = interpolator.get_grid_axis(0).get_cubic_spacing_ratios(ceiling)[2];
    EXPECT_DOUBLE_EQ(result, 1.0);

    result = interpolator.get_grid_axis(1).get_cubic_spacing_ratios(floor)[0];
    EXPECT_DOUBLE_EQ(result, 1.0);
}

TEST_F(CubicImplementationFixture, interpolate)
{
    interpolator.set_target(target);
    std::vector<double> result = interpolator.get_results();
    EXPECT_THAT(result, testing::ElementsAre(testing::DoubleEq(4.158), testing::DoubleEq(11.836)));
}

TEST_F(CubicImplementationFixture, get_cubic_spacing_ratios)
{
    // for cubic dimension 0: {6, 10, 15, 20}, multipliers should be:
    // floor: {4/4, 5/9, 5/10}
    // ceiling: {4/9, 5/10, 5/5}
    std::vector<std::vector<double>> expected_results {{4.0 / 4, 5.0 / 9, 5.0 / 10},
                                                       {4.0 / 9, 5.0 / 10, 5.0 / 5}};
    double result;
    for (std::size_t floor_or_ceiling = 0; floor_or_ceiling <= 1; floor_or_ceiling++) {
        for (std::size_t index = 0; index < 3; index++) {
            result = interpolator.get_grid_axis(0).get_cubic_spacing_ratios(floor_or_ceiling)[index];
            EXPECT_DOUBLE_EQ(result, expected_results[floor_or_ceiling][index]);
        }
    }
}

TEST_F(Grid2DImplementationFixture, consolidate_methods)
{
    interpolator.set_target(target);

    std::vector<Method> expected_methods {Method::linear, Method::linear};
    EXPECT_EQ(interpolator.get_current_methods(), expected_methods);

    std::vector<double> out_of_bounds_vector = {12, 3};
    interpolator.set_target(out_of_bounds_vector);
    expected_methods = {Method::linear, Method::constant};
    EXPECT_EQ(interpolator.get_current_methods(), expected_methods);
}

TEST_F(Grid2DImplementationFixture, get_grid_axis)
{
    std::vector<double> returned_vec = interpolator.get_grid_axis(1).get_values();
    EXPECT_THAT(returned_vec, testing::ElementsAre(4, 6));
}

TEST_F(Grid2DImplementationFixture, get_grid_point_data)
{
    std::vector<std::size_t> coords = {0, 1};
    std::vector<double> returned_vec = interpolator.get_grid_point_data(coords);
    EXPECT_THAT(returned_vec, testing::ElementsAre(3, 6));

    coords = {1, 0};
    returned_vec = interpolator.get_grid_point_data(coords);
    EXPECT_THAT(returned_vec, testing::ElementsAre(2, 4));
}

TEST_F(Grid2DImplementationFixture, get_grid_point_data_relative)
{
    std::vector<std::size_t> coords {0, 1};
    std::vector<short> translation {1, 0}; // {1, 1} stays as is
    std::vector<double> expected_vec = interpolator.get_grid_point_data({1, 1});
    EXPECT_EQ(interpolator.get_grid_point_data_relative(coords, translation), expected_vec);

    translation = {1, 1}; // {1, 2} -> {1, 1}
    EXPECT_EQ(interpolator.get_grid_point_data_relative(coords, translation), expected_vec);

    translation = {-1, 0}; // {-1, 1} -> {0, 1}
    expected_vec = interpolator.get_grid_point_data({0, 1});
    EXPECT_EQ(interpolator.get_grid_point_data_relative(coords, translation), expected_vec);

    translation = {3, -2}; // {3, -1} -> {2, 0}
    expected_vec = interpolator.get_grid_point_data({2, 0});
    EXPECT_EQ(interpolator.get_grid_point_data_relative(coords, translation), expected_vec);
}

} // namespace Btwxt
