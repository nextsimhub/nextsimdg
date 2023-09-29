/*!
 * @file XiosDisabled_test.cpp
 * @brief Initialisation Testing for XIOS
 * @date Feb 28, 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 * 
 * Test file intended to validate the behaviour of the Xios class when Xios is
 * disabled in the parameter file. It does this by setting "enable = false" in 
 * Configurator. 
 * 
 * These unit tests should accompany a corresponding system test where nextsim
 * is run using a config file where under "[xios]"" there is a setting "enable 
 * = false"
 * 
 * These tests should not depend on iodef.xml or any setting in the corresponding
 * environment variable XIOS_IODEF_PATH
 * 
 */

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <mpi.h>
#include "include/Xios.hpp"
#include <iostream>

#include "include/Configurator.hpp"

Nextsim::Xios *xios_handler;

int main( int argc, char* argv[] ) {
    // Disable xios in the 'config'
    Nextsim::Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl
            << "enable = false" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Nextsim::Configurator::addStream(std::move(pcstream));

    // (Temporary) MPI Init - Required by XIOS
    // TODO: Remove this MPI Init
    MPI_Init(NULL, NULL);

    // XIOS Class Init --- Initialises server. Unknown error if initialised per test.
    // TODO: Investigate and find workaround
    xios_handler = new Nextsim::Xios;

    int result = Catch::Session().run( argc, argv );

    // global clean-up...
    delete xios_handler;
    MPI_Finalize();
    return result;
}

// Lock down 'if enabled' protecting Xios::Xios
TEST_CASE("XiosNoInitValidation") 
{
    // Validate XIOS settings when config enabled is set to false
    REQUIRE(xios_handler->validateServerConfiguration());

    // Validate initialization method 
    REQUIRE(xios_handler->validateCalendarConfiguration());
}

// Lock down 'if enabled' protecting XIOS::finalise
TEST_CASE("TestXiosNoFinalise") 
{
    xios_handler->finalise();
}
