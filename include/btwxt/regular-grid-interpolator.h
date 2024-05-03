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
    explicit RegularGridInterpolator(
        const std::vector<GridAxis>& grid_axes);

    RegularGridInterpolator(
        const std::vector<GridAxis>& grid_axes,
        const std::vector<GridPointDataSet>& grid_point_data_sets);

    ~RegularGridInterpolator();

    RegularGridInterpolator(const RegularGridInterpolator& source);

    RegularGridInterpolator& operator=(const RegularGridInterpolator& source);

    void set_axis_extrapolation_method(std::size_t axis_index, ExtrapolationMethod method);

    void set_axis_interpolation_method(std::size_t axis_index, InterpolationMethod method);

    void set_axis_extrapolation_limits(std::size_t axis_index,
                                       const std::pair<double, double>& extrapolation_limits);

    // Public getters
    std::size_t get_number_of_dimensions();

    std::size_t get_number_of_grid_points();

    std::size_t get_number_of_grid_point_data_sets();

    const GridPointDataSet& get_grid_point_data_set(std::size_t data_set_index);

    const GridAxis& get_grid_axis(std::size_t axis_index);

    // Public normalization methods
    double normalize_grid_point_data_set_at_target(std::size_t data_set_index, double scalar = 1.0);

    void normalize_grid_point_data_sets_at_target(double scalar = 1.0);

    // Get results
    void set_target(const std::vector<double>& target);

    double get_value_at_target(std::size_t data_set_index);

    std::vector<double> get_values_at_target();

    [[nodiscard]] std::vector<std::size_t> get_neighboring_indices_at_target() const;

    const std::vector<double>& get_target();

    [[nodiscard]] const std::vector<TargetBoundsStatus>& get_target_bounds_status() const;

  private:
    std::unique_ptr<RegularGridInterpolatorImplementation> implementation;
};

} // namespace Btwxt
