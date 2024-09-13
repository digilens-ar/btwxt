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
    result = interpolator.get_axis_cubic_spacing_ratios(0, floor)[0];
    EXPECT_DOUBLE_EQ(result, 1.0);

    result = interpolator.get_axis_cubic_spacing_ratios(0, ceiling)[0];
    EXPECT_DOUBLE_EQ(result, (10. - 6.) / (15.0 - 6.0));

    result = interpolator.get_axis_cubic_spacing_ratios(0, floor)[1];
    EXPECT_DOUBLE_EQ(result, (15. - 10.) / (15.0 - 6.0));

    result = interpolator.get_axis_cubic_spacing_ratios(0, ceiling)[2];
    EXPECT_DOUBLE_EQ(result, 1.0);

    result = interpolator.get_axis_cubic_spacing_ratios(1, floor)[0];
    EXPECT_DOUBLE_EQ(result, 1.0);
}

TEST_F(CubicImplementationFixture, switch_interp_method)
{
    for (auto i = 0u; i < interpolator.get_number_of_grid_axes(); i++) {
        interpolator.set_axis_interpolation_method(i, InterpolationMethod::cubic);
    }
    interpolator.set_target(target);
    std::vector<double> result1 = interpolator.get_results();
    for (auto i = 0u; i < interpolator.get_number_of_grid_axes(); i++) {
        interpolator.set_axis_interpolation_method(i, InterpolationMethod::linear);
    }
    interpolator.set_target(target); // TODO this shouldn't be necesary but if we don't do this then the changes to the interpolation method won't be taken into account
    std::vector<double> result2 = interpolator.get_results();
    EXPECT_NE(result1, result2);
}

TEST_F(CubicImplementationFixture, interpolate)
{
    interpolator.set_target(target);
    std::vector<double> result = interpolator.get_results();
    EXPECT_THAT(result, testing::ElementsAre(testing::DoubleEq(4.158), testing::DoubleEq(11.836)));
}

TEST_F(CubicImplementationFixture, grid_point_interp_coeffs)
{
    static constexpr std::size_t floor = 0;
    static constexpr std::size_t ceiling = 1;

    interpolator.set_target(target);

    double mu = interpolator.get_floor_to_ceiling_fractions()[0];
    std::size_t floor_grid_point_index = interpolator.get_floor_grid_point_coordinates()[0];

    EXPECT_EQ(interpolator.get_interpolation_coefficients()[0][0],
              2 * mu * mu * mu - 3 * mu * mu + 1);
    EXPECT_EQ(interpolator.get_interpolation_coefficients()[0][1], -2 * mu * mu * mu + 3 * mu * mu);

    EXPECT_EQ(interpolator.get_cubic_slope_coefficients()[0][0],
              (mu * mu * mu - 2 * mu * mu + mu) *
                  interpolator.get_axis_cubic_spacing_ratios(0, floor)[floor_grid_point_index]);
    EXPECT_EQ(interpolator.get_cubic_slope_coefficients()[0][1],
              (mu * mu * mu - mu * mu) *
                  interpolator.get_axis_cubic_spacing_ratios(0, ceiling)[floor_grid_point_index]);
}

TEST_F(CubicImplementationFixture, hypercube_weigh_one_vertex)
{
    interpolator.set_axis_interpolation_method(1, InterpolationMethod::cubic);
    interpolator.set_target(target);
    std::vector<Method> methods = interpolator.get_current_methods();

    std::vector<double> mus = interpolator.get_floor_to_ceiling_fractions();
    double mx = mus[0];
    double my = mus[1];
    double c0x = 2 * mx * mx * mx - 3 * mx * mx + 1;
    double c0y = 2 * my * my * my - 3 * my * my + 1;
    // double c1x = -2*mx*mx*mx + 3*mx*mx;
    double c1y = -2 * my * my * my + 3 * my * my;
    double d0x = mx * mx * mx - 2 * mx * mx + mx;
    double d0y = my * my * my - 2 * my * my + my;
    double d1x = mx * mx * mx - mx * mx;
    double d1y = my * my * my - my * my;
    double s1x = 5.0 / 10;
    double s1y = 2.0 / 4;
    double s0x = 5.0 / 9;
    double s0y = 2.0 / 4;

    std::vector<short> this_vertex = {0, 0};
    double weight = interpolator.get_grid_point_weighting_factor(this_vertex);
    double expected_result = c0x * c0y;
    expected_result += -1 * c0x * d1y * s1y;
    expected_result += -1 * d1x * s1x * c0y;
    expected_result += d1x * s1x * d1y * s1y;
    EXPECT_DOUBLE_EQ(weight, expected_result);

    this_vertex = {-1, 1};
    weight = interpolator.get_grid_point_weighting_factor(this_vertex);
    expected_result = -1 * d0x * s0x * c1y;
    expected_result += -1 * d0x * s0x * d0y * s0y;
    EXPECT_DOUBLE_EQ(weight, expected_result);

    this_vertex = {2, 0};
    weight = interpolator.get_grid_point_weighting_factor(this_vertex);
    expected_result = d1x * s1x * c0y;
    expected_result += -1 * d1x * s1x * d1y * s1y;
    EXPECT_DOUBLE_EQ(weight, expected_result);

    this_vertex = {2, 2};
    weight = interpolator.get_grid_point_weighting_factor(this_vertex);
    expected_result = d1x * s1x * d1y * s1y;
    EXPECT_DOUBLE_EQ(weight, expected_result);
}

TEST_F(CubicImplementationFixture, hypercube_calculations)
{
    interpolator.set_axis_interpolation_method(1, InterpolationMethod::cubic);
    interpolator.set_target(target);

    auto result = interpolator.get_results();
    EXPECT_NEAR(result[0], 4.1953, 0.0001);
    EXPECT_NEAR(result[1], 11.9271, 0.0001);
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
            result = interpolator.get_axis_cubic_spacing_ratios(0, floor_or_ceiling)[index];
            EXPECT_DOUBLE_EQ(result, expected_results[floor_or_ceiling][index]);
        }
    }
}

TEST_F(Grid2DImplementationFixture, grid_point_basics)
{
    interpolator.set_target(target);

    std::vector<std::size_t> expected_floor {1, 0};
    EXPECT_EQ(interpolator.get_floor_grid_point_coordinates(), expected_floor);

    std::vector<double> expected_floor_to_ceiling_fractions {0.4, 0.5};
    EXPECT_EQ(interpolator.get_floor_to_ceiling_fractions(), expected_floor_to_ceiling_fractions);
}

TEST_F(Grid2DImplementationFixture, grid_point_out_of_bounds)
{
    std::vector<double> out_of_bounds_vector = {16, 3};
    interpolator.set_target(out_of_bounds_vector);

    std::vector<std::size_t> expected_floor {1, 0};
    EXPECT_EQ(interpolator.get_floor_grid_point_coordinates(), expected_floor);

    std::vector<double> expected_floor_to_ceiling_fractions {1.2, -0.5};
    EXPECT_EQ(interpolator.get_floor_to_ceiling_fractions(), expected_floor_to_ceiling_fractions);
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

TEST_F(Grid2DImplementationFixture, interpolation_coefficients)
{
    interpolator.set_target(target);

    std::vector<double> mu = interpolator.get_floor_to_ceiling_fractions();

    EXPECT_EQ(interpolator.get_interpolation_coefficients()[0][1], mu[0]);
    EXPECT_EQ(interpolator.get_interpolation_coefficients()[1][0], 1 - mu[1]);

    EXPECT_EQ(interpolator.get_cubic_slope_coefficients()[0][0], 0);
    EXPECT_EQ(interpolator.get_cubic_slope_coefficients()[1][1], 0);
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

TEST_F(Grid3DImplementationFixture, hypercube)
{
    auto& hypercube = interpolator.get_hypercube();
    EXPECT_EQ(hypercube.size(), 16u);
}

TEST_F(Grid3DImplementationFixture, test_hypercube)
{
    auto& hypercube = interpolator.get_hypercube();
    EXPECT_EQ(hypercube.size(), 2u * 4u * 2u);
    EXPECT_THAT(hypercube[0], testing::ElementsAre(0, -1, 0));
    EXPECT_THAT(hypercube[2], testing::ElementsAre(0, 0, 0));
    EXPECT_THAT(hypercube[12], testing::ElementsAre(1, 1, 0));
    EXPECT_THAT(hypercube[15], testing::ElementsAre(1, 2, 1));
}

TEST_F(Grid3DImplementationFixture, make_linear_hypercube)
{
    interpolator.set_axis_interpolation_method(1, InterpolationMethod::linear);
    auto& hypercube = interpolator.get_hypercube();
    EXPECT_EQ(hypercube.size(), 8u);
    EXPECT_THAT(hypercube[0], testing::ElementsAre(0, 0, 0));
    EXPECT_THAT(hypercube[2], testing::ElementsAre(0, 1, 0));
    EXPECT_THAT(hypercube[5], testing::ElementsAre(1, 0, 1));
}

} // namespace Btwxt
