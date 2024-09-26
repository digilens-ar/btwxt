/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Vendor
#include <gtest/gtest.h>

// btwxt
#include <btwxt/btwxt.h>

namespace Btwxt {


inline std::vector<GridAxis> construct_grid_axes(const std::vector<std::vector<double>>& grid_axis_vectors)
{
    std::vector<GridAxis> grid_axes;
    grid_axes.reserve(grid_axis_vectors.size());
    for (const auto& axis : grid_axis_vectors) {
        grid_axes.emplace_back(axis);
    }
    return grid_axes;
}

inline std::vector<std::vector<double>> construct_grid_point_data_sets(const std::vector<std::vector<double>>& grid_point_data_vectors)
{
    std::vector<std::vector<double>> grid_point_data_sets;
    grid_point_data_sets.reserve(grid_point_data_vectors.size());
    for (const auto& grid_point_data_set : grid_point_data_vectors) {
        grid_point_data_sets.emplace_back(
            grid_point_data_set);
    }
    return grid_point_data_sets;
}

class GridImplementationFixture : public testing::Test {
  public:
    std::vector<double> target;
    RegularGridInterpolator interpolator;

    GridImplementationFixture(std::vector<GridAxis> const& gridd, std::vector<std::vector<double>> const& datasetss, InterpolationMethod intMethod):
        interpolator(gridd, construct_grid_point_data_sets(datasetss), intMethod)
    {}

};

class Grid2DImplementationFixture : public GridImplementationFixture {
  protected:
    Grid2DImplementationFixture():
        GridImplementationFixture({GridAxis({0, 10, 15}), GridAxis({4, 6})},  {{6,
                      3,   // 0
                      2,
                      8,   // 10
                      4,
                      2},  // 15
                     {12,
                      6,   // 0
                      4,
                      16,  // 10
                      8,
                      4}}, // 15
            InterpolationMethod::linear)
    {
        target = {12, 5};
    }
};

class CubicImplementationFixture : public GridImplementationFixture {
  protected:
    CubicImplementationFixture():
        GridImplementationFixture({GridAxis({6, 10, 15, 20}), GridAxis({2, 4, 6, 8})},
            {{4,
              3,
              1.5,
              1, // 6
              5,
              4,
              2,
              1, // 10
              8,
              6,
              3,
              2, // 15
              10,
              8,
              4,
              2}, // 20

             {12,
              10,
              4,
              4, // 6
              16,
              12,
              6,
              4, // 10
              20,
              16,
              8,
              4, // 15
              25,
              20,
              10,
              5}}, // 20
            InterpolationMethod::cubic)
    {
        target = {12, 4.5};
    }
};


class Grid3DImplementationFixture : public GridImplementationFixture {
  protected:
    Grid3DImplementationFixture():
        GridImplementationFixture({GridAxis({-15, 0.2, 105}), GridAxis({0, 10, 15}), GridAxis({4, 6})},
           {{6, 3, 2, 8, 4, 2, 3, 6, 13, 2, 0, 15, 3, 6, 13, 2, 0, 15}}, InterpolationMethod::linear)
    {
        target = {26.9, 12, 5};
    }
};

} // namespace Btwxt