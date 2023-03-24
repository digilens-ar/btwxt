/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

// Standard
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <memory>

// btwxt
#include <btwxt/btwxt.h>

// testing
#include "public-fixtures.h"

namespace Btwxt {

TEST(GridAxis, free_check_sorted) {
  std::vector<std::pair<std::vector<double>, bool>> my_vecs = {{{1, 3, 5, 7, 9}, true},
                                                               {{1, 3, 5, 17, 9}, false},
                                                               {{9, 7, 5, 3, 1}, false},
                                                               {{1, 3, 3, 7, 9}, false},
                                                               {{9}, true}};
  bool is_sorted;
  for (const auto &pair : my_vecs) {
    is_sorted = Btwxt::free_check_sorted(pair.first);
    EXPECT_EQ(is_sorted, pair.second);
  }
}

TEST(GridAxis, sorting) {
  std::vector<double> grid_vector = {0, 5, 7, 17, 15};
  auto logger = std::make_shared<BtwxtContextCourierr>();
  EXPECT_THROW(GridAxis my_grid_axis = GridAxis(grid_vector, logger), BtwxtException);
  grid_vector = {0, 5, 7, 10, 15};
  EXPECT_NO_THROW(GridAxis my_grid_axis = GridAxis(grid_vector, logger););
}

TEST(GridAxis, calculate_cubic_spacing_ratios) {
  static constexpr std::size_t floor = 0;
  static constexpr std::size_t ceiling = 1;
  std::vector<double> grid_vector{6, 10, 15, 20, 22};

  GridAxis test_gridaxis(grid_vector, std::make_shared<BtwxtContextCourierr>(), Method::constant,
                         Method::cubic, {-DBL_MAX, DBL_MAX});
  EXPECT_THAT(test_gridaxis.get_cubic_spacing_ratios(floor),
              testing::ElementsAre(1, 5.0 / 9, 0.5, 2.0 / 7));
  EXPECT_THAT(test_gridaxis.get_cubic_spacing_ratios(ceiling),
              testing::ElementsAre(4.0 / 9, 0.5, 5.0 / 7, 1));
}

TEST(GridAxis, bad_limits) {
  auto logger = std::make_shared<BtwxtContextCourierr>();
  GridAxis my_grid_axis({0, 5, 7, 11, 12, 15}, logger);
  std::pair<double, double> extrap_limits{4, 17};
  std::string ExpectedOut = "  NOTE: The lower extrapolation limit (4) is within the set of "
                            "axis values. Setting to smallest axis value (0).\n";
  EXPECT_STDOUT(my_grid_axis.set_extrapolation_limits(extrap_limits);, ExpectedOut);
  EXPECT_EQ(my_grid_axis.get_extrapolation_limits().first, 0);

  extrap_limits = {-2, 12};
  ExpectedOut = "  NOTE: The upper extrapolation limit (12) is within the set of axis values. "
                "Setting to largest axis value (15).\n";
  EXPECT_STDOUT(my_grid_axis.set_extrapolation_limits(extrap_limits);, ExpectedOut);
  EXPECT_EQ(my_grid_axis.get_extrapolation_limits().second, 15);
}
} // namespace Btwxt