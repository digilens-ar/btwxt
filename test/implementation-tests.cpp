/* Copyright (c) 2018 Big Ladder Software LLC. All rights reserved.
 * See the LICENSE file for additional terms and conditions. */

// Standard
#include <chrono>

// vendor
#include <fmt/format.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

// btwxt
#include "fixtures/implementation-fixtures.h"

namespace Btwxt {

TEST_F(CubicImplementationFixture, interpolate)
{
    auto result = interpolator.solve(target);
    EXPECT_THAT(result, testing::ElementsAre(testing::DoubleEq(4.1953125000000018), testing::DoubleEq(11.927125)));
}

TEST_F(Grid2DImplementationFixture, get_grid_axis)
{
    std::vector<double> returned_vec = interpolator.get_grid_axis(1).get_values();
    EXPECT_THAT(returned_vec, testing::ElementsAre(4, 6));
}


} // namespace Btwxt
