/*!
 * @file Logged.hpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_LOGGED_HPP
#define SRC_INCLUDE_LOGGED_HPP

#include <map>
#include <string>
#include <filesystem>

namespace Nextsim {

//! A class to provide general logging facilities.
class Logged {
public:
    //! Static function that configures the logger.
    enum class level {
        ALL,
        TRACE,
        DEBUG,
        INFO,
        NOTICE,
        WARNING,
        ERROR,
        CRITICAL,
        ALERT,
        EMERGENCY,
        NONE
    };

    enum {
        MINIMUM_LOG_LEVEL_KEY,
        FILE_NAME_PATTERN_KEY,
        CONSOLE_LOG_LEVEL_KEY,
    };
    static const std::map<std::string, level> levelNames;
    /*!
     * @brief Initialize the logger by configuring boost::log
     *
     * @param outputDir (Optional) The output directory for storing log files
     */
    static void configure(const std::filesystem::path& outputDir);

    /*!
     * @brief Logs a message at the given log level, or default to level::NOTICE.
     *
     * @param message The message to be logged.
     * @param lvl The level at which to log the message. Defaults to
     *            level::NOTICE if not provided.
     */
    static void log(const std::string& message, const level lvl = level::NOTICE);
    /*!
     * @brief Logs a message at level::TRACE, intended for tracing code execution.
     *
     * @param message The message to be logged.
     */
    static void trace(const std::string& message) { log(message, level::TRACE); };
    /*!
     * @brief Logs a message at level::DEBUG, intended for code debugging.
     *
     * @param message The message to be logged.
     */
    static void debug(const std::string& message) { log(message, level::DEBUG); };
    /*!
     * @brief Logs a message at level::INFO, intended for informational
     * messages that would not normally be shown.
     *
     * @param message The message to be logged.
     */
    static void info(const std::string& message) { log(message, level::INFO); };
    /*!
     * @brief Logs a message at level::NOTICE, intended for messages that would
     * appear during normal execution.
     *
     * @param message The message to be logged.
     */
    static void notice(const std::string& message) { log(message, level::NOTICE); };
    /*!
     * @brief Logs a message at level::WARNING, intended for abnormal
     * conditions that do not affect the continuing execution of the model.
     *
     * @param message The message to be logged.
     */
    static void warning(const std::string& message) { log(message, level::WARNING); };
    /*!
     * @brief Logs a message at level::ERROR, intended for when the model
     * reaches an unrecoverable state.
     *
     * @param message The message to be logged.
     */
    static void error(const std::string& message) { log(message, level::ERROR); };
    /*!
     * @brief Logs a message at level::CRITICAL, intended for critical
     * situations, including detection of problems with the hardware or
     * software environment.
     *
     * @param message The message to be logged.
     */
    static void critical(const std::string& message) { log(message, level::CRITICAL); };
    /*!
     * @brief Logs a message at level::ALERT. Probably not needed for a
     * geophysical model.
     *
     * @param message The message to be logged.
     */
    static void alert(const std::string& message) { log(message, level::ALERT); };
    /*!
     * @brief Logs a message at level::EMERGENCY. Probably not needed for a
     * geophysical model.
     *
     * @param message The message to be logged.
     */
    static void emergency(const std::string& message) { log(message, level::EMERGENCY); };

protected:
    Logged() = default;
    // TODO: Add implementation to actually do some logging
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_LOGGED_HPP */
