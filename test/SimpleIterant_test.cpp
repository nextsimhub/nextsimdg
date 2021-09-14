/*!
 * @file SimpleIterant_test.cpp
 * @date 13 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "SimpleIterant.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

namespace Nextsim {
TEST_CASE( "SimpleIterant testing", "[SimpleIterant]" )
{
    SimpleIterant simps;

    int nSteps = 5;

    Iterator::TimePoint startTime = Iterator::Clock::now();
    Iterator::Duration dt = Iterator::Duration(1);
    Iterator::TimePoint stopTime = startTime + nSteps * dt;

    simps.start(startTime);
    for (int i = 0; i < nSteps; i++) {
        simps.iterate(dt);
    }
    simps.stop(stopTime);
}
}
