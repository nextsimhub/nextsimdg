/*!
 * @file XiosBuild_test.cpp
 * @brief Initialisation Testing for XIOS
 * @date 8th March 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 * 
 * Test file for verifying calls to XIOS class made when XIOS is not
 * linked does not impact operation of Nextsim-DG.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "include/Xios.hpp"

// TODO: Before PR. Consider if this is the best way of running this test. 
// For completeness and diagnosis it could be beneficial however, it may be
// more maintanable to build the XiosInit tests under a different name without
// linking against XIOS.

TEST_CASE("Full XIOS method sweep") 
{
    // Configure the XIOS server
    Nextsim::Xios *xios;
    xios = new Nextsim::Xios;
    xios->initialise();

    // Validate
    REQUIRE(xios->validateConfiguration());
    // Delete XIOS object
    delete xios;   
}
