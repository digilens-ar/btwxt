/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

// Standard
#include <iostream>
#include <numeric>

// btwxt
#include "btwxt.h"
#include "error.h"

namespace Btwxt {

RegularGridInterpolator::RegularGridInterpolator() = default;

RegularGridInterpolator::RegularGridInterpolator(GriddedData &grid_data_in)
    : grid_data(grid_data_in), grid_point(grid_data) {}

RegularGridInterpolator::RegularGridInterpolator(const std::vector<std::vector<double>> &grid,
                                                 const std::vector<std::vector<double>> &values)
    : grid_data(grid, values), grid_point(grid_data) {}

RegularGridInterpolator::RegularGridInterpolator(const RegularGridInterpolator &source) {
  *this = source;
}

std::size_t RegularGridInterpolator::add_value_table(std::vector<double> &value_vector) {
  try {
    return grid_data.add_value_table(value_vector);
  } catch (BtwxtErr &e) {
    log_message(MsgLevel::MSG_ERR, e.what());
    return 0;
  }
}

double RegularGridInterpolator::get_value_at_target(std::vector<double> target,
                                                    std::size_t table_index) {
  set_new_target(target);
  auto [results, err] = grid_point.get_results();
  if (err.has_value()) {
    log_message(MsgLevel::MSG_WARN, err.value());
  }
  return results[table_index];
}

double RegularGridInterpolator::get_value_at_target(std::size_t table_index) {
  auto [results, err] = grid_point.get_results();
  if (err.has_value()) {
    log_message(MsgLevel::MSG_WARN, err.value());
  }
  return results[table_index];
}

std::vector<double>
RegularGridInterpolator::get_values_at_target(const std::vector<double> &target) {
  set_new_target(target);
  auto [results, err] = grid_point.get_results();
  if (err.has_value()) {
    log_message(MsgLevel::MSG_WARN, err.value());
  }
  return results;
}

std::vector<double> RegularGridInterpolator::get_values_at_target() {
  auto [results, err] = grid_point.get_results();
  if (err.has_value()) {
    log_message(MsgLevel::MSG_WARN, err.value());
  }
  return results;
}

double RegularGridInterpolator::normalize_values_at_target(std::size_t table_index,
                                                           const std::vector<double> &target,
                                                           const double scalar) {
  set_new_target(target);
  return normalize_values_at_target(table_index, scalar);
}

double RegularGridInterpolator::normalize_values_at_target(std::size_t table_index,
                                                           const double scalar) {
  try {
    return grid_point.normalize_grid_values_at_target(table_index, scalar);
  } catch (BtwxtWarn &w) {
    log_message(MsgLevel::MSG_WARN, w.what());
    return scalar; // TODO:
  }
}

void RegularGridInterpolator::normalize_values_at_target(const std::vector<double> &target,
                                                         const double scalar) {
  set_new_target(target);
  normalize_values_at_target(scalar);
}

void RegularGridInterpolator::normalize_values_at_target(const double scalar) {
  try {
    return grid_point.normalize_grid_values_at_target(scalar);
  } catch (BtwxtWarn &w) {
    log_message(MsgLevel::MSG_WARN, w.what());
  } catch (BtwxtErr &e) {
    log_message(MsgLevel::MSG_ERR, e.what());
  }
}

void RegularGridInterpolator::set_new_target(const std::vector<double> &target) {
  try {
    grid_point.set_target(target);
  } catch (BtwxtErr &e) {
    log_message(MsgLevel::MSG_ERR, e.what());
  }
}

std::vector<double> RegularGridInterpolator::get_current_target() {
  auto [target, err] = grid_point.get_current_target();
  if (err.has_value()) {
    log_message(MsgLevel::MSG_WARN, err.value());
  }
  return target;
}

void RegularGridInterpolator::clear_current_target() { grid_point = GridPoint(grid_data); }

std::size_t RegularGridInterpolator::get_ndims() { return grid_data.get_ndims(); }

void RegularGridInterpolator::set_axis_extrap_limits(
    const std::size_t &dim, const std::pair<double, double> &extrap_limits) {
  log_message(MsgLevel::MSG_INFO, grid_data.set_axis_extrap_limits(dim, extrap_limits));
}

std::vector<std::vector<short>> &RegularGridInterpolator::get_hypercube() {
  return grid_point.get_hypercube();
}

std::pair<double, double> RegularGridInterpolator::get_axis_limits(int dim) {
  return grid_data.get_extrap_limits(dim);
}

void RegularGridInterpolator::set_logging_callback(BtwxtLoggerFn callback_function,
                                                   void *caller_info) {
  callback_function_ = callback_function;
  caller_context_ = caller_info;
}

void RegularGridInterpolator::log_message(MsgLevel messageType, std::string_view message) {
  if (callback_function_) {
    callback_function_(messageType, message, caller_context_);
  } else if (btwxtCallbackFunction) {
    showMessage(messageType, std::string{message});
  }
}

} // namespace Btwxt
