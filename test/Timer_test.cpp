/*
 * @file Timer_test.cpp
 *
 * @date Oct 27, 2021
 * @author Tim Spain
 */

#include "include/Timer.hpp"


#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <sstream>

TEST_CASE("Test a timer", "[Timer]")
{
    Nextsim::Timer::main.tick("Level 1");
    Nextsim::Timer::main.tick("Level 2a");
    Nextsim::Timer::main.tock("Level 2a");
    Nextsim::Timer::main.tick("Level 2b");
    Nextsim::Timer::main.tock("Level 2b");
    Nextsim::Timer::main.tick("Level 2a");
    Nextsim::Timer::main.tock("Level 2a");
    Nextsim::Timer::main.tock("Level 1");

    std::stringstream sout;
    sout << Nextsim::Timer::main;
    std::cout << sout.str() << std::endl;

    // TODO: Parse the output to actually test something
}
