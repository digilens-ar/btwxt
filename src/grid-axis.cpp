/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <btwxt/grid-axis.h>
#include <spdlog/spdlog.h>

namespace Btwxt {

GridAxis::GridAxis(std::vector<double> values_in,
                   InterpolationMethod int_method)
    :
    values(std::move(values_in))
    , interpolation_method(int_method)
{
    if (values.empty()) {
        throw std::runtime_error("Cannot create grid axis from a zero-length vector.");
    }
    if (!vector_is_valid(values)) {
        throw std::runtime_error("Values are not sorted, or have duplicates.");
    }

    std::fill(cubic_spacing_ratios.begin(), cubic_spacing_ratios.end(), std::vector<double>(values.size() - 1, 1.0));
    if (interpolation_method == InterpolationMethod::cubic) {
        //calculate_cubic_spacing_ratios
        if (values.size() == 1) {
            interpolation_method = InterpolationMethod::linear;
            spdlog::warn("A cubic interpolation method is not valid for grid axis with only one value. "
                         "Interpolation method reset to linear.");
        }
        if (interpolation_method == InterpolationMethod::linear) {
            return;
        }
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

const std::vector<double>& GridAxis::get_cubic_spacing_ratios(size_t floor_or_ceiling) const
{
    return cubic_spacing_ratios[floor_or_ceiling];
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
