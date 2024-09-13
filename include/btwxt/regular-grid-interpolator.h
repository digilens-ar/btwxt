/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Standard
#include <memory>
#include <vector>

// btwxt
#include "grid-axis.h"

namespace Btwxt {

using GridPointDataSet = std::vector<double>; // Data corresponding to all points within a collection of grid axes. Length of data should equal the total number of permutations of grid axes points.

class RegularGridInterpolatorImplementation;

enum class TargetBoundsStatus {
    below_lower_extrapolation_limit,
    extrapolate_low,
    interpolate,
    extrapolate_high,
    above_upper_extrapolation_limit
};

// this will be the public-facing class.
class RegularGridInterpolator {
  public:
    RegularGridInterpolator(
        const std::vector<GridAxis>& grid_axes,
        const std::vector<GridPointDataSet>& grid_point_data_sets);

    ~RegularGridInterpolator();

    RegularGridInterpolator(const RegularGridInterpolator& source);

    RegularGridInterpolator& operator=(const RegularGridInterpolator& source);

    // Public getters
    std::size_t get_number_of_dimensions() const;

    std::size_t get_number_of_grid_points() const;

    std::size_t get_number_of_grid_point_data_sets() const;

    const GridPointDataSet& get_grid_point_data_set(std::size_t data_set_index) const;

    const GridAxis& get_grid_axis(std::size_t axis_index) const;

    // Get results
    void set_target(const std::vector<double>& target);

    double get_value_at_target(std::size_t data_set_index) const;

    std::vector<double> get_values_at_target() const;

    [[nodiscard]] const std::vector<TargetBoundsStatus>& get_target_bounds_status() const;

  private:
    std::unique_ptr<RegularGridInterpolatorImplementation> implementation;
};

} // namespace Btwxt
