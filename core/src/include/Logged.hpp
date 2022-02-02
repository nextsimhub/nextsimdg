/*!
 * @file Logged.hpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_LOGGED_HPP
#define SRC_INCLUDE_LOGGED_HPP

#include <map>
#include <string>

namespace Nextsim {

class Logged {
public:
    //! Static function that configures the logger.
    static void configure();
    enum class level { TRACE, DEBUG, INFO, NOTICE, WARNING, ERROR, CRITICAL, ALERT, EMERGENCY, NONE };

    enum {
        MINIMUM_LOG_LEVEL_KEY,
        FILE_NAME_PATTERN_KEY,
        CONSOLE_LOG_LEVEL_KEY,
    };
    static const std::map<std::string, level> levelNames;

    static void log(const std::string& message, const level lvl);
    static void trace(const std::string& message) { log(message, level::TRACE); };
    static void debug(const std::string& message) { log(message, level::DEBUG); };
    static void info(const std::string& message) { log(message, level::INFO); };
    static void notice(const std::string& message) { log(message, level::NOTICE); };
    static void warning(const std::string& message) { log(message, level::WARNING); };
    static void error(const std::string& message) { log(message, level::ERROR); };
    static void critical(const std::string& message) { log(message, level::CRITICAL); };
    static void alert(const std::string& message) { log(message, level::ALERT); };
    static void emergency(const std::string& message) { log(message, level::EMERGENCY); };

    protected:
    Logged() = default;
    // TODO: Add implementation to actually do some logging
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_LOGGED_HPP */
