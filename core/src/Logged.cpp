/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/Logged.hpp"

#include "include/Configured.hpp"

#include <boost/log/attributes/clock.hpp>
#include <boost/log/common.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <iostream>

using namespace std::literals::string_literals;

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

static constexpr std::array<const char*, static_cast<int>(Logged::level::count)> PRINTED_LEVEL_NAMES
    = { "all", "trace", "debug", "info", "warning", "error", "critical", "fatal", "alert",
          "emergency", "none" };

BOOST_LOG_ATTRIBUTE_KEYWORD(Severity, "Severity", Logged::level)

template <typename CharT, typename TraitsT>
std::basic_ostream<CharT, TraitsT>& operator<<(
    std::basic_ostream<CharT, TraitsT>& stream, const Logged::level& lvl)
{
    stream << PRINTED_LEVEL_NAMES[static_cast<int>(lvl)];
    return stream;
}

boost::log::sources::severity_logger<Logged::level> sl;

// Initialize the logger, that is set up boost::log how we want it
void Logged::configure()
{
    // additional attributes to use during formatting
    boost::log::register_simple_formatter_factory<Logged::level, char>("Severity");
    boost::log::core::get()->add_global_attribute(
        "TimeStamp", boost::log::attributes::local_clock());

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

    // the format string used for both console and file outputs
    constexpr const char* format = "[%TimeStamp%] [%Severity%] %Message%";

    level minimumLogLevel = levelNames.at(
        Configured<Logged>::getConfiguration(keyMap.at(MINIMUM_LOG_LEVEL_KEY), "info"s));
    std::string fileNamePattern
        = Configured<Logged>::getConfiguration(keyMap.at(FILE_NAME_PATTERN_KEY), "nextsim.%T.log"s);

    boost::log::add_file_log(boost::log::keywords::file_name = fileNamePattern,
        // All logs go to file above the minimum level
        boost::log::keywords::filter = (Severity >= minimumLogLevel),
        boost::log::keywords::format = format);

    level consoleLogLevel = levelNames.at(
        Configured<Logged>::getConfiguration(keyMap.at(CONSOLE_LOG_LEVEL_KEY), "error"s));

    boost::log::add_console_log(std::cout,
        boost::log::keywords::filter = (Severity >= consoleLogLevel),
        boost::log::keywords::format = format);
}

void Logged::log(const std::string& message, Logged::level lvl)
{
    BOOST_LOG_SEV(sl, lvl) << message;
}

} /* namespace Nextsim */
