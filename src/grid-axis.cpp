/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

#include <btwxt/grid-axis.h>
#include <spdlog/spdlog.h>

namespace Btwxt {

GridAxis::GridAxis(std::vector<double> values_in,
                   InterpolationMethod interpolation_method,
                   ExtrapolationMethod extrapolation_method,
                   std::pair<double, double> extrapolation_limits)
    :
    values(std::move(values_in))
    , interpolation_method(interpolation_method)
    , extrapolation_method(extrapolation_method)
    , extrapolation_limits(std::move(extrapolation_limits))
    , cubic_spacing_ratios(
          2, std::vector<double>(std::max(static_cast<int>(values.size()) - 1, 0), 1.0))
{
    if (values.empty()) {
        throw std::runtime_error("Cannot create grid axis from a zero-length vector.");
    }
    if (!vector_is_valid(values)) {
        throw std::runtime_error("Values are not sorted, or have duplicates.");
    }
    set_extrapolation_limits(extrapolation_limits);
    set_interpolation_method(interpolation_method);
}

void GridAxis::set_interpolation_method(InterpolationMethod interpolation_method_in)
{
    interpolation_method = interpolation_method_in;
    if (interpolation_method_in == InterpolationMethod::cubic) {
        //calculate_cubic_spacing_ratios
        if (get_length() == 1) {
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

void GridAxis::set_extrapolation_method(ExtrapolationMethod extrapolation_method_in)
{
    switch (extrapolation_method_in) {
    case ExtrapolationMethod::linear: {
        if (get_length() == 1) {
            extrapolation_method = ExtrapolationMethod::constant;
            spdlog::warn(
                "A linear extrapolation method is not valid for grid axis with only one value. "
                "Extrapolation method reset to constant.");
            return;
        }
        break;
    }
    default: {
        break;
    }
    }
    extrapolation_method = extrapolation_method_in;
}

void GridAxis::set_extrapolation_limits(std::pair<double, double> limits)
{
    extrapolation_limits = limits;
    constexpr std::string_view error_format {"{} extrapolation limit ({:.6g}) is within the range "
                                             "of grid axis values [{:.6g}, {:.6g}]."};
    if (extrapolation_limits.first > values[0]) {
        throw std::runtime_error(fmt::format(
            error_format, "Lower", extrapolation_limits.first, values[0], values.back()));
    }
    if (extrapolation_limits.second < values.back()) {
        throw std::runtime_error(fmt::format(
            error_format, "Upper", extrapolation_limits.second, values[0], values.back()));
    }
}

const std::vector<double>& GridAxis::get_cubic_spacing_ratios(bool floor_or_ceiling) const
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
