/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Standard
#include <vector>
#include "grid-axis.h"

namespace Btwxt {


class RegularGridInterpolator {
  public:

    RegularGridInterpolator(const std::vector<GridAxis>& grid_axes,
                                          const std::vector<std::vector<double>>& grid_point_data_sets, InterpolationMethod intMethod=InterpolationMethod::linear);

    [[nodiscard]] std::vector<double> solve(std::vector<double> const& target);

    // Public getters

    [[nodiscard]] std::size_t get_number_of_grid_axes() const
    {
        return grid_axes_.size();
    };

    [[nodiscard]] std::size_t get_number_of_grid_points() const
    {
        return number_of_grid_points_;
    };

    [[nodiscard]] const GridAxis& get_grid_axis(std::size_t axis_index) const
    {
        return grid_axes_[axis_index];
    };

  private:
    // Structured data
    std::vector<GridAxis> grid_axes_;
    std::vector<double> grid_point_data_;
    size_t numDataSets_;
    std::size_t number_of_grid_points_ {0u};
    std::vector<std::size_t> grid_axis_step_size_;   // Used to translate grid point coordinates to
                                                    // indices (size = number_of_grid_axes)

    InterpolationMethod interpolation_method_;
    std::vector<std::array<std::vector<double>, 2>>
        cubic_spacing_ratios_; // Used for cubic interpolation. Outer vector is size 2: 0: spacing
                              // for the floor, 1: spacing for the ceiling. Inner vector is length
                              // of axis values, but the floor vector doesn't use the first entry
                              // and the ceiling doesn't use the last entry.

    // calculated data
    std::vector<std::vector<short>> hypercube; // A minimal set of indices near the target needed to
                                               // perform interpolation calculations.
    std::unordered_map<size_t, std::vector<std::vector<double>>> hypercube_cache;

    // Internal methods
    std::size_t get_grid_point_index_relative(const std::vector<std::size_t>& coordinates,
                                              const std::vector<short>& translation) const;

    std::vector<std::vector<double>> set_hypercube_grid_point_data(
        std::vector<size_t> const& floor_grid_point_coordinates);

    std::vector<double> get_grid_point_data(std::size_t grid_point_index) const;
};

} // namespace Btwxt
