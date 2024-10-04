/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <btwxt/grid-axis.h>
#include <spdlog/spdlog.h>

namespace Btwxt {

GridAxis::GridAxis(std::vector<double> values_in)
    :
    values(std::move(values_in))
{
    if (values.empty()) {
        throw std::runtime_error("Cannot create grid axis from a zero-length vector.");
    }
    if (!vector_is_valid(values)) {
        throw std::runtime_error("Values are not sorted, or have duplicates.");
    }
}


// return an evenly spaced 1-d vector of doubles.
std::vector<double> linspace(double start, double stop, std::size_t number_of_points)
{
    std::vector<double> result(number_of_points);
    double step = (stop - start) / (static_cast<double>(number_of_points) - 1.);
    double value = start;
    for (std::size_t i = 0; i < number_of_points; i++) {
        result[i] = value;
        value += step;
    }
    return result;
}

} // namespace Btwxt
