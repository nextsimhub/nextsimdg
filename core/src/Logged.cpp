/*!
 * @file Logged.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Logged.hpp"

#include <boost/log/common.hpp>
#include <boost/log/utility/setup/file.hpp>

namespace Nextsim {

// Initialize the logger, that is set up boost::log how we want it
void Logged::configure()
{
    boost::log::add_file_log(
            boost::log::keywords::file_name = "nextsim_%Timestep%.log"
            // All logs go to file
    );
}

void Logged::log(const std::string& message, Logged::level lvl)
{
    switch (lvl) {
    case (level::TRACE):
    case (level::DEBUG):
    case (level::INFO):
    case (level::NOTICE):
    case (level::WARNING):
    case (level::ERROR):
    case (level::CRITICAL):
    case (level::ALERT):
        // TODO implement these levels
        break;
    case (level::EMERGENCY):
        break;
    default:
        break;
    }
}

} /* namespace Nextsim */
