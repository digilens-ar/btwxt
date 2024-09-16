/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Standard
#include <memory>
#include <vector>

// btwxt
#include <map>

#include "grid-axis.h"

namespace Btwxt {

using GridPointDataSet = std::vector<double>; // Data corresponding to all points within a collection of grid axes. Length of data should equal the total number of permutations of grid axes points.

class RegularGridInterpolator {
  public:

    RegularGridInterpolator(const std::vector<GridAxis>& grid_axes,
                                          const std::vector<GridPointDataSet>& grid_point_data_sets);

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

    [[nodiscard]] const GridPointDataSet&
    get_grid_point_data_set(std::size_t data_set_index) const
    {
        return grid_point_data_sets_[data_set_index];
    };

    [[nodiscard]] std::size_t get_number_of_grid_point_data_sets() const
    {
        return grid_point_data_sets_.size();
    };

  private:
    // Structured data
    std::vector<GridAxis> grid_axes_;
    std::vector<GridPointDataSet> grid_point_data_sets_;
    std::size_t number_of_grid_points_ {0u};
    std::vector<std::size_t> grid_axis_step_size_;   // Used to translate grid point coordinates to
                                                    // indices (size = number_of_grid_axes)

    // calculated data
    std::vector<std::vector<short>> hypercube; // A minimal set of indices near the target needed to
                                               // perform interpolation calculations.
    std::map<size_t, std::vector<std::vector<double>>> hypercube_cache;

    // Internal methods
    [[nodiscard]] std::size_t
    get_grid_point_index(const std::vector<std::size_t>& coordinates) const;

    std::vector<double> get_grid_point_data_relative(const std::vector<std::size_t>& coordinates,
                                                     const std::vector<short>& translation) const;

    double get_grid_point_weighting_factor(const std::vector<short>& hypercube_indices, std::vector<std::array<double, 4>> const& weighting_factors) const;

    std::size_t get_grid_point_index_relative(const std::vector<std::size_t>& coordinates,
                                              const std::vector<short>& translation) const;

    // for each axis, the fraction the target value
    // is between its floor and ceiling axis values
    std::vector<double> calculate_floor_to_ceiling_fractions(std::vector<double> const& target, std::vector<size_t> const& floor_grid_point_coordinates) const;

    std::vector<std::array<double, 4>> calculate_interpolation_coefficients(std::vector<double> const& floor_to_ceiling_fractions, std::vector<size_t> const& floor_grid_point_coordinates) const;

    std::vector<std::vector<double>> set_hypercube_grid_point_data(
        std::vector<size_t> const& floor_grid_point_coordinates);

    std::vector<double> get_grid_point_data(std::size_t grid_point_index) const;
};

} // namespace Btwxt
