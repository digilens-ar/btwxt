/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef GRID_AXIS_H_
#define GRID_AXIS_H_

#include <algorithm>
#include <array>
#include <vector>


namespace Btwxt {

enum class InterpolationMethod { linear, cubic };
enum class ExtrapolationMethod { constant, linear };

class GridAxis {
    // A single input dimension of the grid
  public:

    explicit GridAxis(
        std::vector<double> values,
        InterpolationMethod interpolation_method = InterpolationMethod::linear,
        ExtrapolationMethod extrapolation_method = ExtrapolationMethod::constant,
        std::pair<double, double> extrapolation_limits = {-DBL_MAX, DBL_MAX});

    // Getters
    [[nodiscard]] const std::vector<double>& get_values() const { return values; }
    [[nodiscard]] std::size_t get_length() const { return values.size(); }
    [[nodiscard]] InterpolationMethod get_interpolation_method() const{ return interpolation_method; }
    [[nodiscard]] ExtrapolationMethod get_extrapolation_method() const { return extrapolation_method; }
    [[nodiscard]] std::pair<double, double> get_extrapolation_limits() const { return extrapolation_limits; }

    [[nodiscard]] const std::vector<double>& get_cubic_spacing_ratios(std::size_t floor_or_ceiling) const;

  private:
    std::vector<double> values;
    InterpolationMethod interpolation_method = InterpolationMethod::linear;
    ExtrapolationMethod extrapolation_method = ExtrapolationMethod::constant;
    std::pair<double, double> extrapolation_limits {-DBL_MAX, DBL_MAX};
    std::array<std::vector<double>, 2>
        cubic_spacing_ratios; // Used for cubic interpolation. Outer vector is size 2: 0: spacing
                              // for the floor, 1: spacing for the ceiling. Inner vector is length
                              // of axis values, but the floor vector doesn't use the first entry
                              // and the ceiling doesn't use the last entry.
};


// free functions

/// @brief Check to see if a vector is valid to be a GridAxis
/// @param vector_in
/// @return true if sorted with no duplicates
inline bool vector_is_valid(const std::vector<double>& vector_in)
{
    if (std::is_sorted(std::begin(vector_in), std::end(vector_in))) {
        // If a vector is sorted, any duplicates will be adjacent to each other
        auto it = std::adjacent_find(vector_in.begin(), vector_in.end());
        if (it == vector_in.end()) {
            return true;
        }
        return false;
    }
    return false;
}

template <typename T>
std::vector<std::vector<T>> cartesian_product(const std::vector<std::vector<T>>& dimension_vectors)
{
    std::vector<std::vector<T>> combinations = {{}};
    for (const auto& list : dimension_vectors) {
        std::vector<std::vector<T>> r;
        for (const auto& x : combinations) {
            for (const auto item : list) {
                r.push_back(x);
                r.back().push_back(item);
            }
        }
        combinations = std::move(r);
    }
    return combinations;
}

inline std::vector<std::vector<double>> cartesian_product(std::vector<GridAxis> const& axes)
{
    std::vector<std::vector<double>> vals;
    vals.reserve(axes.size());
    for (auto const& ax : axes) {
        vals.push_back(ax.get_values());
    }
    return cartesian_product<double>(vals);
}

std::vector<double> linspace(double start, double stop, std::size_t number_of_points);

} // namespace Btwxt

#endif // define GRID_AXIS_H_