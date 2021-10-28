/*
 * @file Timer_test.cpp
 *
 * @date Oct 27, 2021
 * @author Tim Spain
 */

#include "include/Timer.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <chrono>
#include <iostream>
#include <sstream>
#include <thread>

TEST_CASE("Test a timer", "[Timer]")
{

    Nextsim::Timer::main.tick("Level 1");
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    Nextsim::Timer::main.tick("Level 2a");
    std::this_thread::sleep_for(std::chrono::milliseconds(15));
    Nextsim::Timer::main.tock("Level 2a");
    Nextsim::Timer::main.tick("Level 2b");
    std::this_thread::sleep_for(std::chrono::milliseconds(35));
    Nextsim::Timer::main.tick("Level π");
    std::this_thread::sleep_for(std::chrono::milliseconds(11));
    Nextsim::Timer::main.tock("Level π");
    Nextsim::Timer::main.tick("Level e+½");
    std::this_thread::sleep_for(std::chrono::milliseconds(11));
    Nextsim::Timer::main.tock("Level e+½");
    Nextsim::Timer::main.tock("Level 2b");
    Nextsim::Timer::main.tick("Level 2a");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    Nextsim::Timer::main.tick("Level 3");
    std::this_thread::sleep_for(std::chrono::milliseconds(11));
    Nextsim::Timer::main.tock("Level 3");
    Nextsim::Timer::main.tock("Level 2a");
    Nextsim::Timer::main.tock("Level 1");
    Nextsim::Timer::main.tick("Level 1bis");
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    Nextsim::Timer::main.tock("Level 1bis");


    std::stringstream sout;
    sout << Nextsim::Timer::main;
    std::cout << sout.str() << std::endl;

    // TODO: Parse the output to actually test something
}
