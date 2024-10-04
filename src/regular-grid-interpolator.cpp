/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <unordered_map>
#include <cassert>
#include <btwxt/regular-grid-interpolator.hpp>

#include <array>
#include <format>
#include <ranges>
#include <stdexcept>

namespace {
    std::vector<double> derasterData(std::vector<std::vector<double>> const& in)
    {
        std::vector<double> out;
        out.reserve(in.size() * in.at(0).size());

        for (size_t i=0; i<in.at(0).size(); i++) {
            for (auto const& dset : in) {
                out.push_back(dset[i]);
            }
        }
        return out;
    }
}

namespace Btwxt {


copiable_mutex_member::copiable_mutex_member(copiable_mutex_member const& other):
    std::mutex()
{}

copiable_mutex_member& copiable_mutex_member::operator=(copiable_mutex_member const& other)
{
    return *this;
}

RegularGridInterpolator::RegularGridInterpolator(
    std::vector<GridAxis> grid_axes,
    size_t numDataSets,
    std::vector<double> grid_point_data_buffer,
    InterpolationMethod intMethod)
    :
    grid_axes_(std::move(grid_axes))
    , numDataSets_(numDataSets)
    , grid_axis_step_size_(grid_axes_.size()),
    interpolation_method_(intMethod),
    grid_point_data_(std::move(grid_point_data_buffer))
{
    if (interpolation_method_ == InterpolationMethod::cubic) {
        for (auto const& ax : grid_axes_)
        {
            auto& cubic_spacing_ratios = cubic_spacing_ratios_.emplace_back(std::array<std::vector<double>, 2>());
            auto const& values = ax.get_values();
            if (values.size() == 1)
                continue; // A cubic interpolation method is not valid for grid axis with only one value. "
            std::fill(cubic_spacing_ratios.begin(), cubic_spacing_ratios.end(), std::vector<double>(values.size() - 1, 1.0));
            //calculate_cubic_spacing_ratios
            static constexpr std::size_t floor = 0;
            static constexpr std::size_t ceiling = 1;
            for (std::size_t i = 0; i < values.size() - 1; i++) {
                double center_spacing = values[i + 1] - values[i];
                if (i != 0) {
                    cubic_spacing_ratios[floor][i] = center_spacing / (values[i + 1] - values[i - 1]);
                }
                if (i + 2 != values.size()) {
                    cubic_spacing_ratios[ceiling][i] = center_spacing / (values[i + 2] - values[i]);
                }
            }
        }
    }

    // set axis sizes and calculate number of grid points
    number_of_grid_points_ = 1;
    for (std::size_t axis_index = grid_axes_.size(); axis_index-- > 0;) {
        const std::size_t length =
            grid_axes_[axis_index].get_values().size(); // length > 0 ensured by GridAxis constructor
        grid_axis_step_size_[axis_index] = number_of_grid_points_;
        number_of_grid_points_ *= length;
    }

    assert(grid_point_data_.size() == number_of_grid_points_ * numDataSets_);

    hypercube = {{}};
    for (auto const& ax : grid_axes_) {
        bool isCubic = interpolation_method_ == InterpolationMethod::cubic;
        if (ax.get_values().size() == 1)
            isCubic = false; // an axis with size 1 must be treated with linear interpolation
        std::vector<std::vector<short>> r;
        for (const auto& x : hypercube) {
            for (const auto item : (isCubic ?  std::initializer_list<short> {-1, 0, 1, 2} :  std::initializer_list<short> {0, 1})) {
                r.push_back(x);
                r.back().push_back(item);
            }
        }
        hypercube = std::move(r);
    }

    hypercube_cache.resize(number_of_grid_points_,  std::pair {copiable_mutex_member{}, std::vector<size_t>(hypercube.size(), std::numeric_limits<size_t>::max())}); // Use size_t max as an indication of uninitialized cache
}

RegularGridInterpolator::RegularGridInterpolator(std::vector<GridAxis> grid_axes, const std::vector<std::vector<double>>& grid_point_data_sets_, InterpolationMethod intMethod):
    RegularGridInterpolator(std::move(grid_axes), grid_point_data_sets_.size(), derasterData(grid_point_data_sets_), intMethod)
{
    for (const auto& grid_point_data_set : grid_point_data_sets_) {
        if (grid_point_data_set.size() != number_of_grid_points_) {
            throw std::runtime_error(std::format(
            "GridPointDataSet: Size ({}) does not match number of grid points ({}).",
            grid_point_data_set.size(),
            number_of_grid_points_));
        }
    }
}

namespace
{
    std::pmr::vector<std::array<double, 4>> calculate_interpolation_coefficients(
        std::pmr::vector<double> const& floor_to_ceiling_fractions,
        std::pmr::vector<size_t> const& floor_grid_point_coordinates,
        std::vector<GridAxis> const& grid_axes,
        std::vector<std::array<std::vector<double>, 2>> const& cubicSpacingRatios,
        bool cubicInterpolation,
        std::pmr::memory_resource* buffer)
    {
        static constexpr std::size_t floor = 0;
        static constexpr std::size_t ceiling = 1;
        std::pmr::vector<std::array<double, 4>> weighting_factors {std::pmr::polymorphic_allocator<std::array<double, 4>>(buffer)};
        weighting_factors.reserve(grid_axes.size());
        for (std::size_t axis_index = 0; axis_index < grid_axes.size(); axis_index++) {
            double mu = floor_to_ceiling_fractions[axis_index];
            std::array<double, 2> interpolation_coefficients;
            std::array<double, 2> cubic_slope_coefficients;
            const bool isCubic = cubicInterpolation && grid_axes[axis_index].get_values().size() > 1;
            if (isCubic) {
                interpolation_coefficients[floor] = 2 * mu * mu * mu - 3 * mu * mu + 1;
                interpolation_coefficients[ceiling] = -2 * mu * mu * mu + 3 * mu * mu;
                cubic_slope_coefficients[floor] =
                    (mu * mu * mu - 2 * mu * mu + mu) * cubicSpacingRatios[axis_index][floor][floor_grid_point_coordinates[axis_index]];
                cubic_slope_coefficients[ceiling] =
                    (mu * mu * mu - mu * mu) * cubicSpacingRatios[axis_index][ceiling][floor_grid_point_coordinates[axis_index]];
            }
            else {
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

    double get_grid_point_weighting_factor(const std::vector<short>& hypercube_indices, std::pmr::vector<std::array<double, 4>> const& weighting_factors)
    {
        assert(hypercube_indices.size() == weighting_factors.size());
        double weighting_factor = 1.0;
        for (std::size_t axis_index = 0; axis_index < weighting_factors.size(); axis_index++) {
            weighting_factor *= weighting_factors[axis_index][hypercube_indices[axis_index] + 1];
        }
        return weighting_factor;
    }

    // for each axis, the fraction the target value
    // is between its floor and ceiling axis values
    std::pmr::vector<double> calculate_floor_to_ceiling_fractions(
        std::span<const double> target, 
        std::pmr::vector<size_t> const& floor_grid_point_coordinates, 
        std::vector<GridAxis> const& grid_axes,
        std::pmr::memory_resource* buffer)
    {
        assert(target.size() == floor_grid_point_coordinates.size());
        assert(target.size() == grid_axes.size());

        auto compute_fraction = [](double x, double start, double end) -> double
        {
            // how far along an edge is the target?
            return (x - start) / (end - start);
        };

        std::pmr::vector<double> out {std::pmr::polymorphic_allocator<double>(buffer)};
        out.reserve(target.size());
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
}

std::pmr::vector<double> RegularGridInterpolator::solve(std::span<double> target_in, std::pmr::memory_resource* rsrc)
{
    std::array<std::byte, 4096> stackBuffer;
    std::pmr::monotonic_buffer_resource buff_(stackBuffer.data(), stackBuffer.size());
    std::pmr::unsynchronized_pool_resource pool_(&buff_);

    //set_target
    assert(target_in.size() == grid_axes_.size());
    std::pmr::vector<std::size_t> floor_grid_point_coordinates(grid_axes_.size(), 0, std::pmr::polymorphic_allocator<size_t>(&pool_)); // coordinates of the grid point <= target
    for (std::size_t axis_index = 0; axis_index < grid_axes_.size(); axis_index += 1) {
        const auto& axis_values = grid_axes_[axis_index].get_values();
        const int length = static_cast<int>(axis_values.size());
        target_in[axis_index] = std::clamp(target_in[axis_index], axis_values[0], axis_values.back()); // This effectively accomplishes constant extrapolation Forcing the target value to be treated as if it is in bounds
        if (target_in[axis_index] == axis_values.back()) [[unlikely]] {
            floor_grid_point_coordinates[axis_index] =
                std::max(length - 2, 0); // length-2 because that's the left side of the (length-2, length-1) edge.
        }
        else {
            auto upper = std::upper_bound(axis_values.begin(), axis_values.end(), target_in[axis_index]); //TODO this is unnecesearily slow for uniformly spaced values. However, profiling shows that this is only taking ~1% of function time, so probably not worth fixing.
            floor_grid_point_coordinates[axis_index] = upper - axis_values.begin() - 1;
        }
    }
     
    std::pmr::vector<double> floor_to_ceiling_fractions = calculate_floor_to_ceiling_fractions(target_in, floor_grid_point_coordinates, grid_axes_, &pool_);
    auto weighting_factors = calculate_interpolation_coefficients(floor_to_ceiling_fractions, floor_grid_point_coordinates, grid_axes_, cubic_spacing_ratios_, interpolation_method_ == InterpolationMethod::cubic, &pool_);
    auto hypercube_grid_point_data = get_hypercube_grid_data_indices(floor_grid_point_coordinates, &pool_);
    // get results
    std::pmr::vector<double> results(numDataSets_, 0, rsrc);
    for (std::size_t hypercube_index = 0; hypercube_index < hypercube.size(); ++hypercube_index) {
        const double hypercube_weight = get_grid_point_weighting_factor(hypercube[hypercube_index], weighting_factors);
        const size_t hcDataIdx = hypercube_grid_point_data[hypercube_index] * numDataSets_;
        for (std::size_t data_set_index = 0; data_set_index < numDataSets_; ++data_set_index) {
            results[data_set_index] += grid_point_data_[hcDataIdx + data_set_index] * hypercube_weight;
        }
    }
    return results;
}

namespace
{
    std::size_t get_grid_point_index(
        std::pmr::vector<std::size_t> const& coords, std::vector<size_t> const& grid_axis_step_size)
    {
        assert(coords.size() == grid_axis_step_size.size());
        std::size_t grid_point_index = 0;
        for (std::size_t axis_index = 0; axis_index < coords.size(); ++axis_index) {
            grid_point_index += coords[axis_index] * grid_axis_step_size[axis_index];
        }
        return grid_point_index;
    }
}
// Internal calculation methods

std::size_t RegularGridInterpolator::get_grid_point_index_relative(
    const std::pmr::vector<std::size_t>& coords, const std::vector<short>& translation, std::pmr::memory_resource* rsrc) const
{
    std::pmr::vector<std::size_t> temporary_coordinates(grid_axes_.size(), std::pmr::polymorphic_allocator<size_t>(rsrc));
    for (std::size_t axis_index = 0; axis_index < coords.size(); axis_index++) {
        int new_coord = static_cast<int>(coords[axis_index]) + translation[axis_index];
        const int length = static_cast<int>(grid_axes_[axis_index].get_values().size());
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
    return get_grid_point_index(temporary_coordinates, grid_axis_step_size_);
}

std::vector<size_t> const& RegularGridInterpolator::get_hypercube_grid_data_indices(
    std::pmr::vector<size_t> const& floor_grid_point_coordinates, std::pmr::memory_resource* rsrc)
{

    const size_t floor_grid_point_index = get_grid_point_index(floor_grid_point_coordinates, grid_axis_step_size_); // Index of the floor_grid_point_coordinates (used for hypercube caching)
    {
        auto& [cacheItemMutex, cacheItem] = hypercube_cache[floor_grid_point_index];
        std::scoped_lock lock(cacheItemMutex);
        if (cacheItem[0] != std::numeric_limits<size_t>::max()) {
            return cacheItem;
        }
        size_t idx=0;
        for (const auto& v : hypercube) {
            cacheItem[idx++] = get_grid_point_index_relative(floor_grid_point_coordinates, v, rsrc);
        }
        return cacheItem;
    }
}
} // namespace Btwxt
