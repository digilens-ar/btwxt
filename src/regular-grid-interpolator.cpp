/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */


#include "regular-grid-interpolator-implementation.h"

namespace Btwxt {

std::vector<GridAxis> construct_grid_axes(const std::vector<std::vector<double>>& grid_axis_vectors)
{
    std::vector<GridAxis> grid_axes;
    grid_axes.reserve(grid_axis_vectors.size());
    for (const auto& axis : grid_axis_vectors) {
        grid_axes.emplace_back(axis,
                               InterpolationMethod::linear,
                               ExtrapolationMethod::constant,
                               std::pair<double, double> {-DBL_MAX, DBL_MAX});
    }
    return grid_axes;
}

std::vector<GridPointDataSet> construct_grid_point_data_sets(const std::vector<std::vector<double>>& grid_point_data_vectors)
{
    std::vector<GridPointDataSet> grid_point_data_sets;
    grid_point_data_sets.reserve(grid_point_data_vectors.size());
    for (const auto& grid_point_data_set : grid_point_data_vectors) {
        grid_point_data_sets.emplace_back(
            grid_point_data_set);
    }
    return grid_point_data_sets;
}

// Constructors

RegularGridInterpolator::RegularGridInterpolator(const std::vector<GridAxis>& grid)
    : implementation(
          std::make_unique<RegularGridInterpolatorImplementation>(grid, std::vector<GridPointDataSet>{}))
{
}

RegularGridInterpolator::RegularGridInterpolator(
    const std::vector<GridAxis>& grid_axes,
    const std::vector<GridPointDataSet>& grid_point_data_sets)
    : implementation(std::make_unique<RegularGridInterpolatorImplementation>(
          grid_axes, grid_point_data_sets))
{
}

RegularGridInterpolator::~RegularGridInterpolator() = default;

RegularGridInterpolator::RegularGridInterpolator(const RegularGridInterpolator& source)
{
    *this = source;
    this->implementation =
        source.implementation
            ? std::make_unique<RegularGridInterpolatorImplementation>(*source.implementation)
            : nullptr;
}

RegularGridInterpolator& RegularGridInterpolator::operator=(const RegularGridInterpolator& source)
{
    implementation =
        source.implementation
            ? std::make_unique<RegularGridInterpolatorImplementation>(*(source.implementation))
            : nullptr;
    return *this;
}


void RegularGridInterpolator::set_axis_extrapolation_method(const std::size_t axis_index,
                                                            const ExtrapolationMethod method)
{
    implementation->set_axis_extrapolation_method(axis_index, method);
}

void RegularGridInterpolator::set_axis_interpolation_method(const std::size_t axis_index,
                                                            const InterpolationMethod method)
{
    implementation->set_axis_interpolation_method(axis_index, method);
}

void RegularGridInterpolator::set_axis_extrapolation_limits(
    const std::size_t axis_index, const std::pair<double, double>& extrapolation_limits)
{
    implementation->set_axis_extrapolation_limits(axis_index, extrapolation_limits);
}

std::size_t RegularGridInterpolator::get_number_of_dimensions()
{
    return implementation->get_number_of_grid_axes();
}

std::size_t RegularGridInterpolator::get_number_of_grid_points()
{
    return implementation->get_number_of_grid_points();
}

std::size_t RegularGridInterpolator::get_number_of_grid_point_data_sets()
{
    return implementation->get_number_of_grid_point_data_sets();
}

const GridAxis& RegularGridInterpolator::get_grid_axis(std::size_t axis_index)
{
    return implementation->get_grid_axis(axis_index);
}

const GridPointDataSet& RegularGridInterpolator::get_grid_point_data_set(std::size_t data_set_index)
{
    return implementation->get_grid_point_data_set(data_set_index);
}

double RegularGridInterpolator::normalize_grid_point_data_set_at_target(std::size_t data_set_index,
                                                                        const double scalar)
{
    return implementation->normalize_grid_point_data_set_at_target(data_set_index, scalar);
}

void RegularGridInterpolator::normalize_grid_point_data_sets_at_target(const double scalar)
{
    return implementation->normalize_grid_point_data_sets_at_target(scalar);
}

// Public calculation methods
void RegularGridInterpolator::set_target(const std::vector<double>& target)
{
    implementation->set_target(target);
}

double RegularGridInterpolator::get_value_at_target(std::size_t data_set_index)
{
    return implementation->get_results()[data_set_index];
}

std::vector<double> RegularGridInterpolator::get_values_at_target()
{
    return implementation->get_results();
}

std::vector<std::size_t> RegularGridInterpolator::get_neighboring_indices_at_target() const
{
    return implementation->get_neighboring_indices_at_target();
}

const std::vector<double>& RegularGridInterpolator::get_target()
{
    return implementation->get_target();
}

const std::vector<TargetBoundsStatus>& RegularGridInterpolator::get_target_bounds_status() const
{
    return implementation->get_target_bounds_status();
}

} // namespace Btwxt
