/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Standard
#include <algorithm>
#include <cfloat>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>


namespace Btwxt {

class GridPointDataSet {
    // Data corresponding to all points within a collection of grid axes. Length of data should
    // equal the total number of permutations of grid axes points.
  public:
    // Constructors
    GridPointDataSet() = default;

    explicit GridPointDataSet(std::vector<double> data)
        : data(std::move(data))
    {
    }

    std::vector<double> data;
};

} // namespace Btwxt
