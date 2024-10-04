/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

// Standard
#include <iostream>
#include <memory>

// Vendor
#include <gmock/gmock.h>
#include "gmock/gmock-matchers.h"
#include <gtest/gtest.h>

// btwxt
#include <btwxt/btwxt.h>

// testing
#include "fixtures/public-fixtures.h" // EXPECT_STDOUT

namespace Btwxt {

TEST(GridAxis, vector_is_sorted)
{
    std::vector<std::pair<std::vector<double>, bool>> axes = {{{1, 3, 5, 7, 9}, true},
                                                              {{1, 3, 5, 17, 9}, false},
                                                              {{9, 7, 5, 3, 1}, false},
                                                              {{1, 3, 3, 7, 9}, false},
                                                              {{9}, true}};
    bool is_valid;
    for (const auto& pair : axes) {
        is_valid = Btwxt::vector_is_valid(pair.first);
        EXPECT_EQ(is_valid, pair.second);
    }
}

TEST(GridAxis, sorting)
{
    std::vector<double> grid_vector = {0, 5, 7, 17, 15};
    EXPECT_THROW(GridAxis my_grid_axis = GridAxis(grid_vector), std::runtime_error);
    grid_vector = {0, 5, 7, 10, 15};
    EXPECT_NO_THROW(GridAxis my_grid_axis = GridAxis(grid_vector););
}

} // namespace Btwxt