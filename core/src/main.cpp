/*!
 * @file main.cpp
 * @date 10 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#include <iostream>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OASIS
#include <oasis_c.h>
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
    MPI_Comm modelCommunicator = MPI_COMM_WORLD;
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
#ifdef USE_MPI
#ifdef USE_OASIS
        /* We must call these oasis routines before any MPI communication takes place, to make sure
         * we have the right communicator, i.e. modelCommunictor and not MPI_COMM_WORLD. */
        int compID; // Not actually used. Only useful for debugging
        const std::string compName = "nextsim"; // Not useful for any setups we have in mind
        OASIS_CHECK_ERR(oasis_c_init_comp(&compID, compName.c_str(), OASIS_COUPLED));
        OASIS_CHECK_ERR(oasis_c_get_localcomm(&modelCommunicator));
#endif // USE_OASIS
        Nextsim::Model model(modelCommunicator);
#else
        Nextsim::Model model;
#endif
        // Apply the model configuration
        model.configure();
        // Run the Model
        model.run();
    }
#ifdef USE_MPI
#ifdef USE_OASIS
    OASIS_CHECK_ERR(oasis_c_terminate());
#endif
    MPI_Finalize();
#endif

    return 0;
}
