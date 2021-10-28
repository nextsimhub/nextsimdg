/*
 * @file ScopedTimer_test.cpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "../src/include/ScopedTimer.hpp"

#include "include/Timer.hpp"

#include <chrono>
#include <iostream>
#include <thread>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

namespace Nextsim {

void timeAndSleep()
{
    ScopedTimer tas("time and sleep");
    std::this_thread::sleep_for(std::chrono::milliseconds(45));
}

TEST_CASE("Test the scope dependent timer", "[LocalTimer]")
{
    ScopedTimer::setTimerAddress(&Timer::main);
    ScopedTimer::timer().reset();

    ScopedTimer testScopeTimer("test scope timer");
    {
        ScopedTimer localScopeTimer("local scope timer");
        timeAndSleep();
        localScopeTimer.substitute("second scope timer");
        std::this_thread::sleep_for(std::chrono::milliseconds(35));
   }
    timeAndSleep();

    const int nint = 10;
    for (int i = 0; i < nint; ++i) {
        ScopedTimer loop("loop timer");
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    for (int i = 0; i < nint; ++i) {
        ScopedTimer loop("loop timer 2");
        timeAndSleep();
    }

    testScopeTimer.substitute("replacement timer");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));

    std::cout << ScopedTimer::timer() << std::endl;


}

} /* namespace Nextsim */
