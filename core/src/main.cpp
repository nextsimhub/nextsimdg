/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Kacper Kornet <kk562@cam.ac.uk>
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
#include "include/ModelMPI.hpp"
#include "include/NetcdfMetadataConfiguration.hpp"

int main(int argc, char* argv[])
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif // USE_MPI

    int return_code = 0;

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
        MPI_Comm modelCommunicator;
        OASIS_CHECK_ERR(oasis_c_get_localcomm(&modelCommunicator));
        Nextsim::ModelMPI& modelMPI = Nextsim::ModelMPI::getInstance(modelCommunicator);
#else
        Nextsim::ModelMPI& modelMPI = Nextsim::ModelMPI::getInstance(MPI_COMM_WORLD);
#endif // USE_OASIS
#endif
        Nextsim::Model model;
        // Apply the model configuration
        model.configure();
        // Run the Model
        try {
            model.run();
        } catch (const std::exception& e) {
            return_code = -1;
            Nextsim::Logged::error(e.what());
        }
#ifdef USE_MPI
#ifdef USE_OASIS
        OASIS_CHECK_ERR(oasis_c_terminate());
#endif
#endif
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif

    return return_code;
}
