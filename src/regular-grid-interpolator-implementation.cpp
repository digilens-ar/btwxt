/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <unordered_map>
#include <cassert>
#include "regular-grid-interpolator-implementation.h"

#include <array>
#include <format>
#include <stdexcept>

namespace Btwxt {

RegularGridInterpolatorImplementation::RegularGridInterpolatorImplementation(
    const std::vector<GridAxis>& grid_axes,
    const std::vector<GridPointDataSet>& grid_point_data_sets)
    :
    grid_axes(grid_axes)
    , grid_point_data_sets(grid_point_data_sets)
    , grid_axis_step_size(grid_axes.size())
    , methods(grid_axes.size(), Method::undefined)
{
    // set axis sizes and calculate number of grid points
    number_of_grid_points = 1;
    for (std::size_t axis_index = grid_axes.size(); axis_index-- > 0;) {
        const std::size_t length =
            grid_axes[axis_index].get_values().size(); // length > 0 ensured by GridAxis constructor
        grid_axis_step_size[axis_index] = number_of_grid_points;
        number_of_grid_points *= length;
    }

    // Check grid point data set sizes
    for (const auto& grid_point_data_set : grid_point_data_sets) {
        if (grid_point_data_set.size() != number_of_grid_points) {
            throw std::runtime_error(std::format(
            "GridPointDataSet: Size ({}) does not match number of grid points ({}).",
            grid_point_data_set.size(),
            number_of_grid_points));
        }
    }
}

std::vector<double> RegularGridInterpolatorImplementation::solve(const std::vector<double>& target_in)
{
    //set_target
    assert(target_in.size() == grid_axes.size());
    std::vector<std::size_t> floor_grid_point_coordinates(grid_axes.size(), 0); // coordinates of the grid point <= target
    std::vector<TargetBoundsStatus> target_bounds_status(grid_axes.size(), TargetBoundsStatus::interpolate);
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index += 1) {
        const auto& axis_values = grid_axes[axis_index].get_values();
        const int length = static_cast<int>(axis_values.size());
        if (target_in[axis_index] < axis_values[0]) {
            target_bounds_status[axis_index] = TargetBoundsStatus::extrapolate_low;
            floor_grid_point_coordinates[axis_index] = 0;
        }
        else if (target_in[axis_index] > axis_values.back()) {
            target_bounds_status[axis_index] = TargetBoundsStatus::extrapolate_high;
            floor_grid_point_coordinates[axis_index] =
                std::max(length - 2,
                         0); // length-2 because that's the left side of the (length-2, length-1) edge.
        }
        else if (target_in[axis_index] == axis_values.back()) {
            target_bounds_status[axis_index] = TargetBoundsStatus::interpolate;
            floor_grid_point_coordinates[axis_index] =
                std::max(length - 2,
                         0); // length-2 because that's the left side of the (length-2, length-1) edge.
        }
        else {
            target_bounds_status[axis_index] = TargetBoundsStatus::interpolate;
            auto upper = std::upper_bound(axis_values.begin(), axis_values.end(), target_in[axis_index]);
            floor_grid_point_coordinates[axis_index] = upper - axis_values.begin() - 1;
        }
    }
     
    std::vector<double> floor_to_ceiling_fractions = calculate_floor_to_ceiling_fractions(target_in, floor_grid_point_coordinates);
    consolidate_methods(floor_to_ceiling_fractions, target_bounds_status);
    auto weighting_factors = calculate_interpolation_coefficients(floor_to_ceiling_fractions, floor_grid_point_coordinates);
    auto hypercube_grid_point_data = set_hypercube_grid_point_data(floor_grid_point_coordinates);
    // get results
    std::vector<double> results(grid_point_data_sets.size(), 0);
    for (std::size_t hypercube_index = 0; hypercube_index < hypercube.size(); ++hypercube_index) {
        const double hypercube_weight = get_grid_point_weighting_factor(hypercube[hypercube_index], weighting_factors);
        for (std::size_t data_set_index = 0; data_set_index < grid_point_data_sets.size();
             ++data_set_index) {
            results[data_set_index] += hypercube_grid_point_data[hypercube_index][data_set_index] *
                                       hypercube_weight;
        }
    }
    return results;
}

std::vector<double> RegularGridInterpolatorImplementation::get_grid_point_data(std::size_t grid_point_index)
{
    std::vector<double> temporary_grid_point_data(grid_point_data_sets.size(), 0.); 
    for (std::size_t i = 0; i < grid_point_data_sets.size(); ++i) {
        temporary_grid_point_data[i] = grid_point_data_sets[i][grid_point_index];
    }
    return temporary_grid_point_data;
}

std::vector<double>
RegularGridInterpolatorImplementation::get_grid_point_data(const std::vector<std::size_t>& coords)
{
    return get_grid_point_data(get_grid_point_index(coords));
}

std::vector<double> RegularGridInterpolatorImplementation::get_grid_point_data_relative(
    const std::vector<std::size_t>& coords, const std::vector<short>& translation)
{
    return get_grid_point_data(get_grid_point_index_relative(coords, translation));
}

// Internal getter methods
std::vector<Method> RegularGridInterpolatorImplementation::get_interpolation_methods() const
{
    std::vector<Method> interpolation_methods(grid_axes.size());
    static const std::unordered_map<InterpolationMethod, Method> interpolation_method_map {
        {InterpolationMethod::linear, Method::linear}, {InterpolationMethod::cubic, Method::cubic}};

    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        interpolation_methods[axis_index] =
            interpolation_method_map.at(grid_axes[axis_index].get_interpolation_method());
    }
    return interpolation_methods;
}

std::size_t RegularGridInterpolatorImplementation::get_grid_point_index(
    const std::vector<std::size_t>& coords) const
{
    std::size_t grid_point_index = 0;
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); ++axis_index) {
        grid_point_index += coords[axis_index] * grid_axis_step_size[axis_index];
    }
    return grid_point_index;
}

double RegularGridInterpolatorImplementation::get_grid_point_weighting_factor(
    const std::vector<short>& hypercube_indices, std::vector<std::array<double, 4>> const& weighting_factors)
{
    double weighting_factor = 1.0;
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        weighting_factor *= weighting_factors[axis_index][hypercube_indices[axis_index] + 1];
    }
    return weighting_factor;
}

// Internal calculation methods

std::size_t RegularGridInterpolatorImplementation::get_grid_point_index_relative(
    const std::vector<std::size_t>& coords, const std::vector<short>& translation)
{
    std::vector<std::size_t> temporary_coordinates(grid_axes.size());
    for (std::size_t axis_index = 0; axis_index < coords.size(); axis_index++) {
        int new_coord = static_cast<int>(coords[axis_index]) + translation[axis_index];
        const int length = static_cast<int>(grid_axes[axis_index].get_values().size());
        if (new_coord < 0) {
            temporary_coordinates[axis_index] = 0u;
        }
        else if (new_coord >= length) {
            temporary_coordinates[axis_index] = length - 1u;
        }
        else {
            temporary_coordinates[axis_index] = new_coord;
        }
    }
    return get_grid_point_index(temporary_coordinates);
}

std::vector<double> RegularGridInterpolatorImplementation::calculate_floor_to_ceiling_fractions(std::vector<double> const& target, std::vector<size_t> const& floor_grid_point_coordinates) const
{
    auto compute_fraction = [](double x, double start, double end) -> double
    {
        // how far along an edge is the target?
        return (x - start) / (end - start);
    };

    std::vector<double> out;
    out.reserve(grid_axes.size());
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); ++axis_index) {
        auto& axis_values = grid_axes[axis_index].get_values();
        if (axis_values.size() > 1) {
            auto floor_index = floor_grid_point_coordinates[axis_index];
            out.push_back(compute_fraction(
                target[axis_index], axis_values[floor_index], axis_values[floor_index + 1]));
        }
        else {
            out.push_back(1);
        }
    }
    return out;
}

void RegularGridInterpolatorImplementation::consolidate_methods(std::vector<double> const& floor_to_ceiling_fractions, std::vector<
                                                                TargetBoundsStatus> const& target_bounds_status)
// If out of bounds, extrapolate according to prescription
// If outside of extrapolation limits, send a warning and perform constant extrapolation.
{
    std::vector<Method> previous_methods = methods;
    methods = get_interpolation_methods();

    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        switch (target_bounds_status[axis_index]) {
        case TargetBoundsStatus::extrapolate_low:
        case TargetBoundsStatus::extrapolate_high:
            methods[axis_index] = Method::constant;
            break;
        case TargetBoundsStatus::interpolate:
            break;
        }
    }
    reset_hypercube |=
        !std::equal(previous_methods.begin(), previous_methods.end(), methods.begin());
    if (reset_hypercube) {
        set_hypercube(methods, floor_to_ceiling_fractions);
    }
}

void RegularGridInterpolatorImplementation::set_hypercube(std::vector<Method> methods_in, std::vector<double> const& floor_to_ceiling_fractions)
{
    assert(methods_in.size() == grid_axes.size());
    std::size_t previous_size = hypercube.size();
    std::vector<std::vector<int>> options(grid_axes.size(), {0, 1});
    reset_hypercube = false;

    hypercube_size_hash = 0;
    std::size_t digit = 1;
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        if (floor_to_ceiling_fractions[axis_index] == 0.0) {
            options[axis_index] = {0};
            reset_hypercube = true;
        }
        else if (methods_in[axis_index] == Method::cubic) {
            options[axis_index] = {-1, 0, 1, 2};
        }
        hypercube_size_hash += options[axis_index].size() * digit;
        digit *= 10;
    }
    hypercube = {{}};
    for (const auto& list : options) {
        std::vector<std::vector<short>> r;
        for (const auto& x : hypercube) {
            for (const auto item : list) {
                r.push_back(x);
                r.back().push_back(static_cast<short>(item));
            }
        }
        hypercube = std::move(r);
    }
}

std::vector<std::array<double, 4>> RegularGridInterpolatorImplementation::calculate_interpolation_coefficients(std::vector<double> const& floor_to_ceiling_fractions, std::vector<size_t> const& floor_grid_point_coordinates) const
{
    static constexpr std::size_t floor = 0;
    static constexpr std::size_t ceiling = 1;
    std::vector<std::array<double, 4>> weighting_factors;
    weighting_factors.reserve(grid_axes.size());
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        double mu = floor_to_ceiling_fractions[axis_index];
        std::array<double, 2> interpolation_coefficients;
        std::array<double, 2> cubic_slope_coefficients;
        if (methods[axis_index] == Method::cubic) {
            interpolation_coefficients[floor] = 2 * mu * mu * mu - 3 * mu * mu + 1;
            interpolation_coefficients[ceiling] = -2 * mu * mu * mu + 3 * mu * mu;
            cubic_slope_coefficients[floor] =
                (mu * mu * mu - 2 * mu * mu + mu) *
                grid_axes[axis_index].get_cubic_spacing_ratios(floor)[floor_grid_point_coordinates[axis_index]];
            cubic_slope_coefficients[ceiling] =
                (mu * mu * mu - mu * mu) *
                grid_axes[axis_index].get_cubic_spacing_ratios(ceiling)[floor_grid_point_coordinates[axis_index]];
        }
        else {
            if (methods[axis_index] == Method::constant) {
                mu = mu < 0 ? 0 : 1;
            }
            interpolation_coefficients[floor] = 1 - mu;
            interpolation_coefficients[ceiling] = mu;
            cubic_slope_coefficients[floor] = 0.0;
            cubic_slope_coefficients[ceiling] = 0.0;
        }
        weighting_factors.emplace_back(std::array<double, 4> {
             -cubic_slope_coefficients[floor], // point below floor (-1)
            interpolation_coefficients[floor] - cubic_slope_coefficients[ceiling], // floor (0)
            interpolation_coefficients[ceiling] + cubic_slope_coefficients[floor], // ceiling (1)
            cubic_slope_coefficients[ceiling] // point above ceiling (2)
        });
    }
    return weighting_factors;
}

std::vector<std::vector<double>>
RegularGridInterpolatorImplementation::set_hypercube_grid_point_data(
    std::vector<size_t> const& floor_grid_point_coordinates)
{
    const size_t floor_grid_point_index = get_grid_point_index(floor_grid_point_coordinates); // Index of the floor_grid_point_coordinates (used for hypercube caching)
    if (hypercube_cache.count({floor_grid_point_index, hypercube_size_hash})) {
        return hypercube_cache.at({floor_grid_point_index, hypercube_size_hash});
    }
    std::size_t hypercube_index = 0;
    std::vector<std::vector<double>> hypercube_grid_point_data(hypercube.size(), std::vector<double>(grid_point_data_sets.size()));
    for (const auto& v : hypercube) {
        hypercube_grid_point_data[hypercube_index] =
            get_grid_point_data_relative(floor_grid_point_coordinates, v);
        ++hypercube_index;
    }
    hypercube_cache[{floor_grid_point_index, hypercube_size_hash}] = hypercube_grid_point_data;
    return hypercube_grid_point_data;
}
} // namespace Btwxt