/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

#include <map>
#include <vector>
#include <string>
#include "btwxt/regular-grid-interpolator.h"



namespace Btwxt {

enum class Method { undefined, constant, linear, cubic };

class RegularGridInterpolatorImplementation {
    friend class GridAxis;

  public:

    RegularGridInterpolatorImplementation(const std::vector<GridAxis>& grid_axes,
                                          const std::vector<GridPointDataSet>& grid_point_data_sets);

    [[nodiscard]] std::vector<double> solve(std::vector<double> const& target);

    // Public getters

    [[nodiscard]] std::size_t get_number_of_grid_axes() const
    {
        return grid_axes.size();
    };

    [[nodiscard]] std::size_t get_number_of_grid_points() const
    {
        return number_of_grid_points;
    };

    [[nodiscard]] const GridAxis& get_grid_axis(std::size_t axis_index) const
    {
        check_axis_index(axis_index, "get grid axis");
        return grid_axes[axis_index];
    };

    [[nodiscard]] const GridPointDataSet&
    get_grid_point_data_set(std::size_t data_set_index) const
    {
        check_data_set_index(data_set_index, "get grid point data set");
        return grid_point_data_sets[data_set_index];
    };

    [[nodiscard]] std::size_t get_number_of_grid_point_data_sets() const
    {
        return grid_point_data_sets.size();
    };

    std::vector<double> get_grid_point_data(const std::vector<std::size_t>& coordinates);

    std::vector<double> get_grid_point_data_relative(const std::vector<std::size_t>& coordinates,
                                                     const std::vector<short>& translation);

    [[nodiscard]] std::size_t
    get_grid_point_index(const std::vector<std::size_t>& coordinates) const;

    double get_grid_point_weighting_factor(const std::vector<short>& hypercube_indices, std::vector<std::array<double, 4>> const& weighting_factors);

  private:
    // Structured data
    std::vector<GridAxis> grid_axes;
    std::vector<GridPointDataSet> grid_point_data_sets;
    std::size_t number_of_grid_points {0u};
    std::vector<std::size_t> grid_axis_step_size;   // Used to translate grid point coordinates to
                                                    // indices (size = number_of_grid_axes)

    // calculated data

    std::vector<Method> methods;
    std::vector<std::vector<short>> hypercube; // A minimal set of indices near the target needed to
                                               // perform interpolation calculations.
    bool reset_hypercube {false};

    std::vector<std::vector<double>> hypercube_grid_point_data;
    std::vector<double> hypercube_weights;

    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::vector<double>>> hypercube_cache;

    std::size_t hypercube_size_hash {0u};

    // Internal methods
    std::size_t get_grid_point_index_relative(const std::vector<std::size_t>& coordinates,
                                              const std::vector<short>& translation);

    // for each axis, the fraction the target value
    // is between its floor and ceiling axis values
    std::vector<double> calculate_floor_to_ceiling_fractions(std::vector<double> const& target, std::vector<size_t> const& floor_grid_point_coordinates) const;

    void consolidate_methods(std::vector<double> const& floor_to_ceiling_fractions, std::vector<double> const& target, std::vector<
                             TargetBoundsStatus> const& target_bounds_status);

    std::vector<std::array<double, 4>> calculate_interpolation_coefficients(std::vector<double> const& floor_to_ceiling_fractions, std::vector<size_t> const& floor_grid_point_coordinates) const;

    void set_hypercube(std::vector<Method> methods, std::vector<double> const& floor_to_ceiling_fractions);

    void set_hypercube_grid_point_data(std::vector<size_t> const& floor_grid_point_coordinates);

    std::vector<double> set_results(std::vector<std::array<double, 4>> const& weighting_factors);

    void check_axis_index(std::size_t axis_index, const std::string& action_description) const;

    void check_data_set_index(std::size_t data_set_index,
                              const std::string& action_description) const;

    [[nodiscard]] std::vector<Method> get_interpolation_methods() const;

    std::vector<double> get_grid_point_data(std::size_t grid_point_index);

};

} // namespace Btwxt
