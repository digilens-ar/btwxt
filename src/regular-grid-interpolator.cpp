/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */


#include "regular-grid-interpolator-implementation.h"

namespace Btwxt {


// Constructors
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

std::size_t RegularGridInterpolator::get_number_of_dimensions() const
{
    return implementation->get_number_of_grid_axes();
}

std::size_t RegularGridInterpolator::get_number_of_grid_points() const
{
    return implementation->get_number_of_grid_points();
}

std::size_t RegularGridInterpolator::get_number_of_grid_point_data_sets() const
{
    return implementation->get_number_of_grid_point_data_sets();
}

const GridAxis& RegularGridInterpolator::get_grid_axis(std::size_t axis_index) const
{
    return implementation->get_grid_axis(axis_index);
}
 
const GridPointDataSet& RegularGridInterpolator::get_grid_point_data_set(std::size_t data_set_index) const
{
    return implementation->get_grid_point_data_set(data_set_index);
}

// Public calculation methods
void RegularGridInterpolator::set_target(const std::vector<double>& target)
{
    implementation->set_target(target);
}

double RegularGridInterpolator::get_value_at_target(std::size_t data_set_index) const
{
    return implementation->get_results()[data_set_index];
}

std::vector<double> RegularGridInterpolator::get_values_at_target() const
{
    return implementation->get_results();
}

const std::vector<TargetBoundsStatus>& RegularGridInterpolator::get_target_bounds_status() const
{
    return implementation->get_target_bounds_status();
}

} // namespace Btwxt
