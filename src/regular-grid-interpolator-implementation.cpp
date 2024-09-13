/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <unordered_map>
#include <cassert>
#include "regular-grid-interpolator-implementation.h"
#include <format>
#include <stdexcept>

namespace Btwxt {

RegularGridInterpolatorImplementation::RegularGridInterpolatorImplementation(
    const std::vector<GridAxis>& grid_axes,
    const std::vector<GridPointDataSet>& grid_point_data_sets)
    :
    grid_axes(grid_axes)
    , grid_point_data_sets(grid_point_data_sets)
    , grid_axis_lengths(grid_axes.size())
    , grid_axis_step_size(grid_axes.size())
    , target(grid_axes.size(), 0.)
    , floor_grid_point_coordinates(grid_axes.size(), 0)
    , target_bounds_status(grid_axes.size())
    , methods(grid_axes.size(), Method::undefined)
    , weighting_factors(grid_axes.size(), std::vector<double>(4, 0.))
    , results(grid_point_data_sets.size(), 0.)
    , interpolation_coefficients(grid_axes.size(), std::vector<double>(2, 0.))
    , cubic_slope_coefficients(grid_axes.size(), std::vector<double>(2, 0.))
{
    // set axis sizes and calculate number of grid points
    number_of_grid_points = 1;
    for (std::size_t axis_index = grid_axes.size(); axis_index-- > 0;) {
        std::size_t length =
            grid_axes[axis_index].get_length(); // length > 0 ensured by GridAxis constructor
        grid_axis_lengths[axis_index] = length;
        grid_axis_step_size[axis_index] = number_of_grid_points;
        number_of_grid_points *= length;
    }

    // Check grid point data set sizes
    for (const auto& grid_point_data_set : grid_point_data_sets) {
        check_grid_point_data_set_size(grid_point_data_set);
    }
}

void RegularGridInterpolatorImplementation::set_target(const std::vector<double>& target_in)
{
    if (target_in.size() != grid_axes.size()) {
        throw std::runtime_error(std::format("Target (size={}) and grid (size={}) do not have the same dimensions.",
                        target_in.size(),
                        grid_axes.size()));
    }
    if (target_is_set) {
        if ((target_in == target) && (methods == get_interpolation_methods())) {
            return;
        }
    }
    target = target_in;
    target_is_set = true;
    set_floor_grid_point_coordinates();
    std::vector<double> floor_to_ceiling_fractions = calculate_floor_to_ceiling_fractions();
    consolidate_methods(floor_to_ceiling_fractions);
    calculate_interpolation_coefficients(floor_to_ceiling_fractions);
    set_results();
}

const std::vector<double>& RegularGridInterpolatorImplementation::get_target() const
{
    if (!target_is_set) {
         throw std::runtime_error("The current target was requested, but no target has been set.");
    }
    return target;
}

std::vector<double> RegularGridInterpolatorImplementation::get_results() const
{
    if (grid_point_data_sets.size() == 0u) {
        throw std::runtime_error("There are no grid point data sets. No results returned.");
    }
    if (!target_is_set) {
         throw std::runtime_error("Results were requested, but no target has been set.");
    }
    return results;
}

void RegularGridInterpolatorImplementation::normalize_grid_point_data_sets_at_target(
    const double scalar)
{
    if (!target_is_set) {
         throw std::runtime_error("Cannot normalize grid point data sets. No target has been set.");
    }
    for (std::size_t data_set_index = 0; data_set_index < grid_point_data_sets.size();
         ++data_set_index) {
        normalize_grid_point_data_set(data_set_index, results[data_set_index] * scalar);
    }
    hypercube_cache.clear();
    set_results();
}

double RegularGridInterpolatorImplementation::normalize_grid_point_data_set_at_target(
    std::size_t data_set_index, double scalar)
{
    check_data_set_index(data_set_index, "normalize grid point data set");
    if (!target_is_set) {
        throw std::runtime_error(std::format(
            "GridPointDataSet '{}': Cannot normalize grid point data set. No target has been set.",
            data_set_index));
    }
    // create a scalar which represents the product of the inverted normalization factor and the
    // value in the data set at the independent variable reference value
    double total_scalar = results[data_set_index] * scalar;
    normalize_grid_point_data_set(data_set_index, total_scalar);
    hypercube_cache.clear();
    set_results();

    return total_scalar;
}

void RegularGridInterpolatorImplementation::normalize_grid_point_data_set(
    std::size_t data_set_index, double scalar)
{
    check_data_set_index(data_set_index, "normalize grid point data set");
    auto& data_set = grid_point_data_sets[data_set_index];
    if (scalar == 0.0) {
         throw std::runtime_error(std::format("GridPointDataSet '{}': Attempt to normalize grid point data set by zero.",
                        data_set_index));
    }
    scalar = 1.0 / scalar;
    std::transform(data_set.begin(),
                   data_set.end(),
                   data_set.begin(),
                   [scalar](double x) -> double { return x * scalar; });
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
    const std::vector<short>& hypercube_indices)
{
    double weighting_factor = 1.0;
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        weighting_factor *= weighting_factors[axis_index][hypercube_indices[axis_index] + 1];
    }
    return weighting_factor;
}

std::vector<std::size_t> RegularGridInterpolatorImplementation::get_neighboring_indices_at_target(std::vector<double> const& floor_to_ceiling_fractions) const
{
    if (!target_is_set) {
        throw std::runtime_error("Cannot retrieve neighboring indices. No target has been set.");
    }
    std::vector<std::vector<std::size_t>> axes_neighbor_indices(
        grid_axes.size(),
        std::vector<std::size_t>()); // For each axis, what are the neighboring indices?
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); ++axis_index) {
        auto floor_index = floor_grid_point_coordinates[axis_index];
        if (floor_to_ceiling_fractions[axis_index] < 1.0) {
            axes_neighbor_indices[axis_index].push_back(floor_index);
        }
        if (grid_axis_lengths[axis_index] > 1 && floor_to_ceiling_fractions[axis_index] > 0.0) {
            axes_neighbor_indices[axis_index].push_back(floor_index + 1);
        }
    }
    std::vector<std::vector<std::size_t>> axes_neighbor_coordinates =
        cartesian_product(axes_neighbor_indices);
    std::vector<std::size_t> neighbor_indices;
    neighbor_indices.reserve(axes_neighbor_coordinates.size());
    for (const auto& coordinates : axes_neighbor_coordinates) {
        neighbor_indices.push_back(get_grid_point_index(coordinates));
    }
    return neighbor_indices;
}

// private methods

void RegularGridInterpolatorImplementation::check_grid_point_data_set_size(
    const GridPointDataSet& grid_point_data_set)
{
    if (grid_point_data_set.size() != number_of_grid_points) {
         throw std::runtime_error(std::format(
            "GridPointDataSet: Size ({}) does not match number of grid points ({}).",
            grid_point_data_set.size(),
            number_of_grid_points));
    }
}

void RegularGridInterpolatorImplementation::set_results()
{
    set_hypercube_grid_point_data();
    std::fill(results.begin(), results.end(), 0.0);
    for (std::size_t hypercube_index = 0; hypercube_index < hypercube.size(); ++hypercube_index) {
        hypercube_weights[hypercube_index] =
            get_grid_point_weighting_factor(hypercube[hypercube_index]);
        for (std::size_t data_set_index = 0; data_set_index < grid_point_data_sets.size();
             ++data_set_index) {
            results[data_set_index] += hypercube_grid_point_data[hypercube_index][data_set_index] *
                                       hypercube_weights[hypercube_index];
        }
    }
}

// Internal calculation methods

std::size_t RegularGridInterpolatorImplementation::get_grid_point_index_relative(
    const std::vector<std::size_t>& coords, const std::vector<short>& translation)
{
    std::vector<std::size_t> temporary_coordinates(grid_axes.size());
    int new_coord;
    for (std::size_t axis_index = 0; axis_index < coords.size(); axis_index++) {
        new_coord = static_cast<int>(coords[axis_index]) + translation[axis_index];
        if (new_coord < 0) {
            temporary_coordinates[axis_index] = 0u;
        }
        else if (new_coord >= static_cast<int>(grid_axis_lengths[axis_index])) {
            temporary_coordinates[axis_index] = grid_axis_lengths[axis_index] - 1u;
        }
        else {
            temporary_coordinates[axis_index] = new_coord;
        }
    }
    return get_grid_point_index(temporary_coordinates);
}

void RegularGridInterpolatorImplementation::set_floor_grid_point_coordinates()
{
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index += 1) {
        set_axis_floor_grid_point_index(axis_index);
    }
    floor_grid_point_index = get_grid_point_index(floor_grid_point_coordinates);
}

void RegularGridInterpolatorImplementation::set_axis_floor_grid_point_index(std::size_t axis_index)
{
    const auto& axis_values = grid_axes[axis_index].get_values();
    int length = static_cast<int>(grid_axis_lengths[axis_index]);
    if (target[axis_index] < get_extrapolation_limits(axis_index).first) {
        target_bounds_status[axis_index] = TargetBoundsStatus::below_lower_extrapolation_limit;
        floor_grid_point_coordinates[axis_index] = 0u;
    }
    else if (target[axis_index] > get_extrapolation_limits(axis_index).second) {
        target_bounds_status[axis_index] = TargetBoundsStatus::above_upper_extrapolation_limit;
        floor_grid_point_coordinates[axis_index] =
            std::max(length - 2,
                     0); // length-2 because that's the left side of the (length-2, length-1) edge.
    }
    else if (target[axis_index] < axis_values[0]) {
        target_bounds_status[axis_index] = TargetBoundsStatus::extrapolate_low;
        floor_grid_point_coordinates[axis_index] = 0;
    }
    else if (target[axis_index] > axis_values.back()) {
        target_bounds_status[axis_index] = TargetBoundsStatus::extrapolate_high;
        floor_grid_point_coordinates[axis_index] =
            std::max(length - 2,
                     0); // length-2 because that's the left side of the (length-2, length-1) edge.
    }
    else if (target[axis_index] == axis_values.back()) {
        target_bounds_status[axis_index] = TargetBoundsStatus::interpolate;
        floor_grid_point_coordinates[axis_index] =
            std::max(length - 2,
                     0); // length-2 because that's the left side of the (length-2, length-1) edge.
    }
    else {
        target_bounds_status[axis_index] = TargetBoundsStatus::interpolate;
        auto upper = std::upper_bound(axis_values.begin(), axis_values.end(), target[axis_index]);
        floor_grid_point_coordinates[axis_index] = upper - axis_values.begin() - 1;
    }
}

void RegularGridInterpolatorImplementation::check_axis_index(std::size_t axis_index,
    const std::string& action_description) const
{
    if (axis_index > grid_axes.size() - 1) {
         throw std::runtime_error(std::format(
            "Axis index, {}, does not exist. Unable to {}. Number of grid axes = {}.",
            axis_index,
            action_description,
            grid_axes.size()));
    }
}

void RegularGridInterpolatorImplementation::check_data_set_index(std::size_t data_set_index,
    const std::string& action_description) const
{
    if (data_set_index > grid_point_data_sets.size() - 1) {
        throw std::runtime_error(std::format("Data set index, {}, does not exist. Unable to {}. Number of "
                               "grid point data sets = {}.",
                               data_set_index,
                               action_description,
                               grid_point_data_sets.size()));
    }
}

std::vector<double> RegularGridInterpolatorImplementation::calculate_floor_to_ceiling_fractions() const
{
    auto compute_fraction = [](double x, double start, double end) -> double
    {
        // how far along an edge is the target?
        return (x - start) / (end - start);
    };

    std::vector<double> out;
    out.reserve(grid_axes.size());
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); ++axis_index) {
        if (grid_axis_lengths[axis_index] > 1) {
            auto& axis_values = grid_axes[axis_index].get_values();
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

void RegularGridInterpolatorImplementation::consolidate_methods(std::vector<double> const& floor_to_ceiling_fractions)
// If out of bounds, extrapolate according to prescription
// If outside of extrapolation limits, send a warning and perform constant extrapolation.
{
    std::vector<Method> previous_methods = methods;
    methods = get_interpolation_methods();
    if (target_is_set) {
        // get extrapolation methods
         std::vector<Method> extrapolation_methods(grid_axes.size());
        static const std::unordered_map<ExtrapolationMethod, Method> extrapolation_method_map {
            {ExtrapolationMethod::constant, Method::constant},
            {ExtrapolationMethod::linear, Method::linear}};
        for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
            extrapolation_methods[axis_index] =
                extrapolation_method_map.at(grid_axes[axis_index].get_extrapolation_method());
        }

        constexpr std::string_view error_format {
            "GridAxis '{}': The target ({:.6g}) is {} the extrapolation "
            "limit ({:.6g})."};
        for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
            switch (target_bounds_status[axis_index]) {
            case TargetBoundsStatus::extrapolate_low:
            case TargetBoundsStatus::extrapolate_high:
                methods[axis_index] = extrapolation_methods[axis_index];
                break;
            case TargetBoundsStatus::below_lower_extrapolation_limit:
                 throw std::runtime_error(std::format(error_format,
                                       axis_index,
                                       target[axis_index],
                                       "below",
                                       get_extrapolation_limits(axis_index).first));
                break;
            case TargetBoundsStatus::above_upper_extrapolation_limit:
                 throw std::runtime_error(std::format(error_format,
                                       axis_index,
                                       target[axis_index],
                                       "above",
                                       get_extrapolation_limits(axis_index).second));
                break;
            case TargetBoundsStatus::interpolate:
                break;
            }
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
        if (target_is_set && floor_to_ceiling_fractions[axis_index] == 0.0) {
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
    if (hypercube.size() != previous_size) {
        hypercube_grid_point_data.resize(hypercube.size(),
                                         std::vector<double>(grid_point_data_sets.size()));
        hypercube_weights.resize(hypercube.size());
    }
}

void RegularGridInterpolatorImplementation::calculate_interpolation_coefficients(std::vector<double> const& floor_to_ceiling_fractions)
{
    static constexpr std::size_t floor = 0;
    static constexpr std::size_t ceiling = 1;
    for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
        double mu = floor_to_ceiling_fractions[axis_index];
        if (methods[axis_index] == Method::cubic) {
            interpolation_coefficients[axis_index][floor] = 2 * mu * mu * mu - 3 * mu * mu + 1;
            interpolation_coefficients[axis_index][ceiling] = -2 * mu * mu * mu + 3 * mu * mu;
            cubic_slope_coefficients[axis_index][floor] =
                (mu * mu * mu - 2 * mu * mu + mu) *
                get_axis_cubic_spacing_ratios(axis_index,
                                              floor)[floor_grid_point_coordinates[axis_index]];
            cubic_slope_coefficients[axis_index][ceiling] =
                (mu * mu * mu - mu * mu) *
                get_axis_cubic_spacing_ratios(axis_index,
                                              ceiling)[floor_grid_point_coordinates[axis_index]];
        }
        else {
            if (methods[axis_index] == Method::constant) {
                mu = mu < 0 ? 0 : 1;
            }
            interpolation_coefficients[axis_index][floor] = 1 - mu;
            interpolation_coefficients[axis_index][ceiling] = mu;
            cubic_slope_coefficients[axis_index][floor] = 0.0;
            cubic_slope_coefficients[axis_index][ceiling] = 0.0;
        }
        weighting_factors[axis_index][0] =
            -cubic_slope_coefficients[axis_index][floor]; // point below floor (-1)
        weighting_factors[axis_index][1] =
            interpolation_coefficients[axis_index][floor] -
            cubic_slope_coefficients[axis_index][ceiling]; // floor (0)
        weighting_factors[axis_index][2] =
            interpolation_coefficients[axis_index][ceiling] +
            cubic_slope_coefficients[axis_index][floor]; // ceiling (1)
        weighting_factors[axis_index][3] =
            cubic_slope_coefficients[axis_index][ceiling]; // point above ceiling (2)
    }
}

void RegularGridInterpolatorImplementation::set_hypercube_grid_point_data()
{
    if (hypercube_cache.count({floor_grid_point_index, hypercube_size_hash})) {
        hypercube_grid_point_data =
            hypercube_cache.at({floor_grid_point_index, hypercube_size_hash});
        return;
    }
    std::size_t hypercube_index = 0;
    for (const auto& v : hypercube) {
        hypercube_grid_point_data[hypercube_index] =
            get_grid_point_data_relative(floor_grid_point_coordinates, v);
        ++hypercube_index;
    }
    hypercube_cache[{floor_grid_point_index, hypercube_size_hash}] = hypercube_grid_point_data;
}
} // namespace Btwxt