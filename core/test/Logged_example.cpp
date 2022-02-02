/*!
 * @file Logged_example.cpp
 *
 * @date Feb 1, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include <sstream>

#include "include/Logged.hpp"
#include "include/Configurator.hpp"
#include "include/Configured.hpp"

int main()
{

    Nextsim::Configurator::clear();

    std::stringstream config;

    config << "[Logged]" << std::endl;
    config << "minimum_log_level = emergency" << std::endl;
//    config << "file_name_pattern = loggy_mclogface.log" << std::endl;
    config << "console_log_level = trace" << std::endl;
    Nextsim::Configurator::addSStream(config);

    Nextsim::Logged::configureLogging();

    Nextsim::Logged::trace("A TRACE message");
    Nextsim::Logged::debug("A DEBUG message");
    Nextsim::Logged::info("An INFO message");
    Nextsim::Logged::warning("A WARNING message");
    Nextsim::Logged::error("An ERROR message");
    Nextsim::Logged::critical("A CRITICAL message");
    Nextsim::Logged::alert("An ALERT message");
    Nextsim::Logged::emergency("An EMERGENCY message");

}

