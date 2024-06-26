/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Vendor
#include <gtest/gtest.h>

// btwxt
#include "regular-grid-interpolator-implementation.h"
#include <btwxt/btwxt.h>

namespace Btwxt {

class GridImplementationFixture : public testing::Test {
  public:
    std::vector<std::vector<double>> grid;
    std::vector<std::vector<double>> data_sets;
    std::vector<double> target;
    RegularGridInterpolatorImplementation interpolator;

    GridImplementationFixture(std::vector<std::vector<double>> const& gridd, std::vector<std::vector<double>> const& datasetss):
        grid(gridd),
        data_sets(datasetss),
        interpolator(construct_grid_axes(grid), construct_grid_point_data_sets(data_sets))
    {}

};

class Grid2DImplementationFixture : public GridImplementationFixture {
  protected:
    Grid2DImplementationFixture():
        GridImplementationFixture({{0, 10, 15}, {4, 6}},  {{6,
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
                      4}}) // 15
    {
        target = {12, 5};
        interpolator.set_axis_extrapolation_method(0, ExtrapolationMethod::linear);
    }
};

class CubicImplementationFixture : public GridImplementationFixture {
  protected:
    CubicImplementationFixture():
        GridImplementationFixture({{6, 10, 15, 20}, {2, 4, 6, 8}},
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
              5}}) // 20
    {
        target = {12, 4.5};
        interpolator.set_axis_interpolation_method(0, InterpolationMethod::cubic);
    }
};


class Grid3DImplementationFixture : public GridImplementationFixture {
  protected:
    Grid3DImplementationFixture():
        GridImplementationFixture({{-15, 0.2, 105}, {0, 10, 15}, {4, 6}},
           {{6, 3, 2, 8, 4, 2, 3, 6, 13, 2, 0, 15, 3, 6, 13, 2, 0, 15}})
    {
        target = {26.9, 12, 5};
        interpolator.set_axis_interpolation_method(0, InterpolationMethod::linear);
        interpolator.set_axis_interpolation_method(1, InterpolationMethod::cubic);
        interpolator.set_axis_interpolation_method(2, InterpolationMethod::linear);
    }
};

} // namespace Btwxt