/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#ifndef GRID_AXIS_H_
#define GRID_AXIS_H_

#include <algorithm>
#include <array>
#include <vector>


namespace Btwxt {

enum class InterpolationMethod { linear, cubic };

class GridAxis {
    // A single input dimension of the grid
  public:
    explicit GridAxis(std::vector<double> values);
    // Getters
    [[nodiscard]] const std::vector<double>& get_values() const { return values; }
  private:
    std::vector<double> values;
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