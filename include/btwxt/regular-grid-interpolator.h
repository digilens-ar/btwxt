/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#pragma once

// Standard
#include <memory_resource>
#include <vector>
#include "grid-axis.h"

namespace Btwxt {

// A mutex that creates a new mutex when you try to copy it. This just allows having a unique mutex class member for each copy of a class without having to declare a copy constructor for that class.
class copiable_mutex_member : public std::mutex
{
public:
	copiable_mutex_member() = default;

	copiable_mutex_member([[maybe_unused]] copiable_mutex_member const& other);

	copiable_mutex_member& operator=([[maybe_unused]] copiable_mutex_member const& other);
};

class RegularGridInterpolator {
  public:

    // grid_point_data_buffer should have its data stored such that the dataset axis is the "fast axis", grid_axes.back() is the 2nd fastest axis, and grid_axes[0] is the "slow axis"
    RegularGridInterpolator(
        std::vector<GridAxis> const& grid_axes,
        size_t numDataSets,
        std::vector<double> grid_point_data_buffer, 
        InterpolationMethod intMethod=InterpolationMethod::linear);

    //This constructor will have some overhead to convert grid_point_data_sets to the form used in the first constructor
    // each entry to grid_point_data_sets should have its data stored such that grid_axes.back() is the "fast axis" and grid_axes[0] is the "slow axis"
    RegularGridInterpolator(
        std::vector<GridAxis> const& grid_axes,
        std::vector<std::vector<double>> const& grid_point_data_sets, 
        InterpolationMethod intMethod=InterpolationMethod::linear);

    // If a value of target is outside the interpolation range, it will modified to lie within range.
    [[nodiscard]] std::pmr::vector<double> solve(std::vector<double>& target, std::pmr::memory_resource* rsrc=std::pmr::get_default_resource());

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
    std::vector<std::vector<short>> hypercube; // The minimal set of indices relative to the target's nearest grid point needed to
                                               // perform interpolation calculations.
    std::vector<std::pair<copiable_mutex_member, std::vector<size_t>>> hypercube_cache; // stores the grid point data indices for each element of the hypercube for a given floor index. Frankly this doesn't seem necessary.

    // Internal methods
    std::size_t get_grid_point_index_relative(const std::pmr::vector<std::size_t>& coordinates,
                                              const std::vector<short>& translation, std::pmr::memory_resource* rsrc) const;

    std::vector<size_t> const& get_hypercube_grid_data_indices(
        std::pmr::vector<size_t> const& floor_grid_point_coordinates, std::pmr::memory_resource* rsrc);
};

} // namespace Btwxt
