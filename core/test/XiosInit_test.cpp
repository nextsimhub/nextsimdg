/*!
 * @file XiosInit_test.cpp
 * @brief Initialisation Testing for XIOS
 * @date Feb 28, 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 * 
 * Test file intended to validate the initialization of the server process of
 * the XIOS class, including syncronising the state of the Nextsim grid. This
 * file serves as a set of unit tests, the system tests are elsewhere.
 */


#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "include/Xios.hpp"

TEST_CASE("Init XIOS, check default validation") 
{
    // Configure the XIOS server
    Nextsim::Xios *xios;
    xios = new Nextsim::Xios::Xios(NULL, NULL);
    xios->initialise();

    // Validate initialization
    REQUIRE(xios->validateServerConfiguration());

    // Validate initialization
    REQUIRE(xios->validateCalendarConfiguration());
    // Delete XIOS object
    delete xios;   
}

// TEST_CASE("Init XIOS, check config validation") 
// {
//     // Configure the XIOS server
//     Nextsim::Xios *xios;
//     xios = new Nextsim::Xios::Xios(NULL, NULL);
//     xios->initialise();

//     std::string xios_timestep = xios->getCalendarTimestep();
//     std::string config_timestep = "P0-0T0:10:0";

//     REQUIRE(xios_timestep == config_timestep);

//     std::string xios_startDate = xios->getCalendarStart();
//     std::string config_startDate = "2010-01-01T10:00:00Z";

//     Nextsim::Xios::convertXiosDateStringToIsoDate(xios_startDate);
//     // Delete XIOS object
//     delete xios;   
// }