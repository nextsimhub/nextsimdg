/*!
 * @file main.cpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#include <iostream>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "include/CommandLineParser.hpp"
#include "include/ConfigurationHelpPrinter.hpp"
#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/Model.hpp"
#include "include/NetcdfMetadataConfiguration.hpp"

int main(int argc, char* argv[])
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif // USE_MPI

    // Pass the command line to Configurator to handle
    Nextsim::Configurator::setCommandLine(argc, argv);
    // Extract any config files defined on the command line
    Nextsim::CommandLineParser cmdLine(argc, argv);
    // Pass the config file names to Configurator
    Nextsim::Configurator::addFiles(cmdLine.getConfigFileNames());
    try {
        // Get the configuration stored in the restart file
        Nextsim::NetcdfMetadataConfiguration ncdfMC;
        Nextsim::Configurator::setAdditionalConfiguration(&ncdfMC);
        Nextsim::Configurator::getAdditionalConfiguration(Nextsim::Model::restartOptionName);
    } catch (const std::exception& e) {
        // Do nothing. If there is no additional configuration to be parse, ignore it.
    }
    // Parse the configuration to load those that are explicitly configured
    Nextsim::ConfiguredModule::parseConfigurator();

    if (!cmdLine.configHelp().empty()) {
        Nextsim::Model::HelpMap map;
        Nextsim::Model::getHelpRecursive(map, true);
        Nextsim::ConfigurationHelpPrinter::setOutput(
            Nextsim::ConfigurationHelpPrinter::Output::ANSI);
        Nextsim::ConfigurationHelpPrinter::print(std::cout, map, cmdLine.configHelp());
    } else {
        // Construct the Model
        Nextsim::Model model;
        // Apply the model configuration
        model.configure();
        // Run the Model
        model.run();
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif

    std::cout << "SUCCESS\n";
    return 0;
}
