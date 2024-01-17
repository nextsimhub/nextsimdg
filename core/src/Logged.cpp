/*!
 * @file Logged.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Logged.hpp"

#include "include/Configured.hpp"

#include <boost/log/common.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <iostream>

namespace Nextsim {

const std::map<std::string, Logged::level> Logged::levelNames = {
    { "all", level::ALL },
    { "All", level::ALL },
    { "ALL", level::ALL },
    { "trace", level::TRACE },
    { "TRACE", level::TRACE },
    { "debug", level::DEBUG_LVL },
    { "DEBUG", level::DEBUG_LVL },
    { "info", level::INFO },
    { "INFO", level::INFO },
    { "warning", level::WARNING },
    { "WARNING", level::WARNING },
    { "error", level::ERROR },
    { "ERROR", level::ERROR },
    { "critical", level::CRITICAL },
    { "CRITICAL", level::CRITICAL },
    { "fatal", level::CRITICAL },
    { "FATAL", level::CRITICAL },
    { "alert", level::ALERT },
    { "ALERT", level::ALERT },
    { "emergency", level::EMERGENCY },
    { "EMERGENCY", level::EMERGENCY },
    { "none", level::NONE },
    { "None", level::NONE },
    { "NONE", level::NONE },
};

BOOST_LOG_ATTRIBUTE_KEYWORD(Severity, "Severity", Logged::level)

boost::log::sources::severity_logger<Logged::level> sl;

// Initialize the logger, that is set up boost::log how we want it
void Logged::configure()
{
    /*
     * The enum and keyMap that would usually be inherited from Configured are
     * defined here so that a class can be derived from Logged
     * without having to be derived from Configured and Boost::program_options.
     */
    enum {
        MINIMUM_LOG_LEVEL_KEY,
        FILE_NAME_PATTERN_KEY,
        CONSOLE_LOG_LEVEL_KEY,
    };

    const std::map<int, std::string> keyMap = {
        { MINIMUM_LOG_LEVEL_KEY, "Logged.minimum_log_level" },
        { FILE_NAME_PATTERN_KEY, "Logged.file_name_pattern" },
        { CONSOLE_LOG_LEVEL_KEY, "Logged.console_log_level" },
    };

    level minimumLogLevel = levelNames.at(Configured<Logged>::getConfiguration(
        keyMap.at(MINIMUM_LOG_LEVEL_KEY), std::string("info")));
    std::string fileNamePattern = Configured<Logged>::getConfiguration(
        keyMap.at(FILE_NAME_PATTERN_KEY), std::string("nextsim.%T.log"));
    boost::log::add_file_log(boost::log::keywords::file_name = fileNamePattern,
        // All logs go to file above the minimum level
        boost::log::keywords::filter = (Severity >= minimumLogLevel));

    level consoleLogLevel = levelNames.at(Configured<Logged>::getConfiguration(
        keyMap.at(CONSOLE_LOG_LEVEL_KEY), std::string("none")));
    boost::log::add_console_log(
        std::cout, boost::log::keywords::filter = (Severity >= consoleLogLevel));
}

void Logged::log(const std::string& message, Logged::level lvl)
{
    BOOST_LOG_SEV(sl, lvl) << message;
}

} /* namespace Nextsim */
