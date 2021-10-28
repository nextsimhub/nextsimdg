/*
 * @file LocalTimer_test.cpp
 *
 * @date Oct 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/LocalTimer.hpp"
#include "include/Timer.hpp"

#include <chrono>
#include <iostream>
#include <thread>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

namespace Nextsim {

void timeAndSleep()
{
    LocalTimer tas("time and sleep");
    std::this_thread::sleep_for(std::chrono::milliseconds(45));
}

TEST_CASE("Test the scope dependent timer", "[LocalTimer]")
{
    LocalTimer::setTimerAddress(&Timer::main);
    LocalTimer::timer().reset();

    LocalTimer testScopeTimer("test scope timer");
    {
        LocalTimer localScopeTimer("local scope timer");
        std::this_thread::sleep_for(std::chrono::milliseconds(35));
        timeAndSleep();
   }
    timeAndSleep();

    const int nint = 10;
    for (int i = 0; i < nint; ++i) {
        LocalTimer loop("loop timer");
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }

    for (int i = 0; i < nint; ++i) {
        LocalTimer loop("loop timer 2");
        timeAndSleep();
    }

    std::cout << LocalTimer::timer() << std::endl;


}

} /* namespace Nextsim */
