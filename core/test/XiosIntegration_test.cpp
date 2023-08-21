/*!
 * @file XiosIntegration_test.cpp
 * @brief Integration/System Testing for XIOS
 * @date 8th March 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 * 
 * Testing if XIOS interfaces correctly with the rest of Nextsim-DG. Enables
 * XIOS via the config file. Running tests for each type of grid and comparing
 * to the expected value
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "include/Xios.hpp"

// TODO: This will need to be added to the overall system test (like custom output dir)
// core/src/ will be an inappropriate place to run this test due to incorrect componentisation

TEST_CASE("XIOS enabled via config file for Dev Grid (config_simple_example)") 
{
    // TODO: Make this a common method (perhaps in custom output dir change).
    // We want to use a DYNAMIC_SECTION over config files here too!

    // // Setup the command line arguments using this test class
    // Nextsim::ArgV argv({ "nextsimdg", "--config-file", targetConfigFilename });

    // // Setup the configurator
    // Nextsim::Configurator::setCommandLine(argv.argc(), argv());
    // Nextsim::CommandLineParser cmdLine(argv.argc(), argv());
    // Nextsim::Configurator::addFiles(cmdLine.getConfigFileNames());
    // Nextsim::ConfiguredModule::parseConfigurator();

    // // Setup the nextsim model
    // Nextsim::Model model;
    // model.configure();
    // model.run();

}